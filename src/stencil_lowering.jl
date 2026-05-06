"""
    stencil_lowering.jl — lower stencil-form discretization rules to ESS
    `replacement` form.

The ESD rule catalog ships rules in two forms:

- **Replacement form** — a single `replacement` ExpressionNode AST that ESS
  `rule_engine` (RFC §7) consumes directly via `parse_rule`.
- **Stencil form** — an array of `{selector, coeff}` entries plus a
  `combine` op. The array is a compact, grid-family-aware representation
  the ESS rule engine does NOT consume.

ESS retired its stencil-aware MMS evaluator (esm-4t5) in service of the
single-pathway directive: every simulation must flow through
`discretize → ArrayOp → eval`. `discretize` only sees rules in
replacement form. Stencil-only catalog entries (`upwind_1st`,
`lax_friedrichs_flux`, …) therefore need a structural lowering before
they can drive the canonical pipeline.

This module owns that lowering. It is a pure AST → AST transform — it
performs no numerical evaluation, runs no simulation, and applies no
per-rule-shape logic. The lowered output is the same JSON-shaped
`Dict{String,Any}` the input was, with a `replacement` field inserted.

Public surface:

- [`lower_stencil_to_replacement`](@ref) — single-rule lowerer.

Selector kinds currently supported: `cartesian`, `arakawa`. Other
families (`vertical`, `latlon`, `cubed_sphere`, `indirect`, `reduction`)
raise `ArgumentError` until their lowering rows are added; each new
kind composes here as a separate dispatch branch.
"""

"""
    lower_stencil_to_replacement(rule::AbstractDict) -> Dict{String,Any}

Return a copy of `rule` with a `replacement` field inserted. Idempotent:
if `rule` already carries `replacement`, the input is returned via
shallow copy (no recomputation).

Stencil-form lowering produces

    combine_op(coeff_i * index(operand, axis_args_i...)
               for entry i in stencil)

where:

- `combine_op` is the rule's top-level `combine` field (defaults to `"+"`).
- `operand` is the first `\$`-prefixed positional argument of `applies_to.args`.

The shape of `axis_args_i` depends on the stencil family:

- **cartesian** (`kind="cartesian"`): single-axis. `axis_pat` comes from
  `applies_to.dim` (must itself be a `\$`-prefixed pattern variable);
  each entry's `selector.axis` must equal `axis_pat`. Lowered form is
  `index(operand, axis_pat + offset_i)`.

- **arakawa** (`kind="arakawa"`): multi-axis. The set of unique
  `selector.axis` pattern variables across the stencil defines the
  positional dimension order in the lowered `index` (sorted
  lexicographically — e.g. `\$x` before `\$y`). Each entry contributes
  `index(operand, \$a + off_a, \$b + off_b, …)` where the entry's `axis`
  slot carries `\$axis + offset` and other slots carry the bare axis
  pattern variable. The selector's `stagger` (`cell_center`, `face_x`,
  `face_y`, `vertex`) is NOT encoded in the lowered AST — it stays
  recoverable via the rule dict's preserved `stencil` field, which the
  binding-side runner uses (per SELECTOR_KINDS.md decision #16) to pick
  the matching component array on a C-grid (`face_x` ↔ x-component,
  `face_y` ↔ y-component). `applies_to.dim` is NOT required for arakawa
  rules — axis pattern variables are intrinsic to the stencil entries.

A stencil with a single entry lowers to the bare term — no `combine_op`
wrapper is emitted. An offset of `0` lowers to `axis_pat` (the trivial
`+ 0` is omitted) at that axis slot.

Errors:

- `ArgumentError` if `rule` has neither `stencil` nor `replacement`.
- `ArgumentError` if `applies_to` or `applies_to.args` is missing or
  malformed.
- `ArgumentError` if entries declare conflicting selector kinds.
- `ArgumentError` if any selector kind is not yet supported.
- `ArgumentError` if a cartesian rule is missing `applies_to.dim` or
  has an entry whose axis disagrees with `applies_to.dim`.
- `ArgumentError` if an arakawa entry's `axis` is not a `\$`-prefixed
  pattern variable, or its `stagger` is not one of the four allowed
  symbols.
- `ArgumentError` if `combine` is not one of `"+"`, `"*"`, `"min"`, `"max"`.
"""
function lower_stencil_to_replacement(rule::AbstractDict)::Dict{String, Any}
    out = Dict{String, Any}(String(k) => v for (k, v) in rule)

    if haskey(out, "replacement")
        return out
    end

    haskey(out, "stencil") || throw(
        ArgumentError(
            "lower_stencil_to_replacement: rule has neither 'replacement' nor 'stencil'",
        ),
    )

    applies_to = get(out, "applies_to", nothing)
    applies_to isa AbstractDict || throw(
        ArgumentError(
            "lower_stencil_to_replacement: rule missing 'applies_to' object",
        ),
    )

    operand_var = _stencil_operand_pattern_var(applies_to)

    stencil = out["stencil"]
    stencil isa AbstractVector || throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'stencil' must be an array",
        ),
    )
    isempty(stencil) && throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'stencil' must contain at least one entry",
        ),
    )

    combine_op = String(get(out, "combine", "+"))
    combine_op in ("+", "*", "min", "max") || throw(
        ArgumentError(
            "lower_stencil_to_replacement: unsupported combine '$combine_op' " *
                "(expected one of '+', '*', 'min', 'max')",
        ),
    )

    family = _stencil_family(stencil)

    terms = if family == "cartesian"
        axis_var = _stencil_axis_pattern_var(applies_to)
        Any[_lower_cartesian_entry_unpack(entry, operand_var, axis_var, i)
            for (i, entry) in enumerate(stencil)]
    elseif family == "arakawa"
        axis_vars = _arakawa_collect_axis_vars(stencil)
        Any[_lower_arakawa_entry_unpack(entry, operand_var, axis_vars, i)
            for (i, entry) in enumerate(stencil)]
    else
        throw(
            ArgumentError(
                "lower_stencil_to_replacement: stencil selector kind " *
                    "'$family' is not yet supported (currently supported: " *
                    "'cartesian', 'arakawa')",
            ),
        )
    end

    replacement = length(terms) == 1 ? terms[1] :
                  Dict{String, Any}("op" => combine_op, "args" => terms)

    out["replacement"] = replacement
    return out
end

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

function _stencil_operand_pattern_var(applies_to::AbstractDict)::String
    args = get(applies_to, "args", nothing)
    args isa AbstractVector || throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'applies_to.args' must be an array",
        ),
    )
    for a in args
        if a isa AbstractString && startswith(String(a), "\$")
            return String(a)
        end
    end
    throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'applies_to.args' must contain at " *
                "least one '\$'-prefixed pattern variable (operand)",
        ),
    )
end

function _stencil_axis_pattern_var(applies_to::AbstractDict)::String
    dim = get(applies_to, "dim", nothing)
    dim isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'applies_to.dim' must be a string",
        ),
    )
    s = String(dim)
    startswith(s, "\$") || throw(
        ArgumentError(
            "lower_stencil_to_replacement: 'applies_to.dim' must be a " *
                "'\$'-prefixed pattern variable (got '$s')",
        ),
    )
    return s
end

"""
    _stencil_family(stencil) -> String

Inspect every entry's `selector.kind` and return the (single) family name
shared across the stencil. All entries must agree; mixed-kind stencils
are not supported.
"""
function _stencil_family(stencil::AbstractVector)::String
    family = ""
    for (i, entry) in enumerate(stencil)
        entry isa AbstractDict || throw(
            ArgumentError(
                "lower_stencil_to_replacement: stencil entry $i must be an object",
            ),
        )
        selector = get(entry, "selector", nothing)
        selector isa AbstractDict || throw(
            ArgumentError(
                "lower_stencil_to_replacement: stencil entry $i missing 'selector' object",
            ),
        )
        kind = String(get(selector, "kind", ""))
        if isempty(family)
            family = kind
        elseif kind != family
            throw(
                ArgumentError(
                    "lower_stencil_to_replacement: stencil mixes selector kinds " *
                        "('$family' and '$kind' at entry $i); a single rule must " *
                        "use one selector family throughout",
                ),
            )
        end
    end
    return family
end

function _entry_selector_and_coeff(entry, idx::Int)
    entry isa AbstractDict || throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx must be an object",
        ),
    )
    selector = get(entry, "selector", nothing)
    selector isa AbstractDict || throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx missing 'selector' object",
        ),
    )
    coeff = get(entry, "coeff", nothing)
    coeff === nothing && throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx missing 'coeff'",
        ),
    )
    return selector, coeff
end

function _lower_cartesian_entry_unpack(entry, operand_var::String, axis_var::String, idx::Int)
    selector, coeff = _entry_selector_and_coeff(entry, idx)
    return _lower_cartesian_entry(selector, coeff, operand_var, axis_var, idx)
end

function _lower_arakawa_entry_unpack(entry, operand_var::String,
        axis_vars::Vector{String}, idx::Int)
    selector, coeff = _entry_selector_and_coeff(entry, idx)
    return _lower_arakawa_entry(selector, coeff, operand_var, axis_vars, idx)
end

function _lower_cartesian_entry(
        selector::AbstractDict,
        coeff,
        operand_var::String,
        axis_var::String,
        idx::Int,
    )
    sel_axis_raw = get(selector, "axis", nothing)
    sel_axis_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: cartesian entry $idx 'selector.axis' " *
                "must be a string",
        ),
    )
    sel_axis = String(sel_axis_raw)
    if sel_axis != axis_var
        throw(
            ArgumentError(
                "lower_stencil_to_replacement: cartesian entry $idx 'selector.axis' " *
                    "= '$sel_axis' disagrees with 'applies_to.dim' = '$axis_var'",
            ),
        )
    end

    offset_raw = get(selector, "offset", nothing)
    (offset_raw isa Integer && !(offset_raw isa Bool)) || throw(
        ArgumentError(
            "lower_stencil_to_replacement: cartesian entry $idx 'selector.offset' " *
                "must be an integer",
        ),
    )
    offset = Int(offset_raw)

    index_arg = if offset == 0
        axis_var
    else
        Dict{String, Any}("op" => "+", "args" => Any[axis_var, offset])
    end

    index_node = Dict{String, Any}(
        "op" => "index", "args" => Any[operand_var, index_arg],
    )

    return Dict{String, Any}("op" => "*", "args" => Any[coeff, index_node])
end

# -----------------------------------------------------------------------------
# Arakawa lowering
# -----------------------------------------------------------------------------

const _ARAKAWA_STAGGERS = ("cell_center", "face_x", "face_y", "vertex")

"""
    _arakawa_collect_axis_vars(stencil) -> Vector{String}

Walk every entry's `selector.axis` and return the sorted set of unique
`\$`-prefixed axis pattern variables. The sort order pins the positional
dimension order of the lowered `index` op (so e.g. `[\$x, \$y]` means
`index(operand, \$x_arg, \$y_arg)`).

Each entry's `axis` MUST be a `\$`-prefixed pattern variable; literal axis
names (e.g. `"lon"` / `"lat"` / `"xi"`) belong to other selector families
and are rejected by the per-entry validator.
"""
function _arakawa_collect_axis_vars(stencil::AbstractVector)::Vector{String}
    seen = Set{String}()
    for (i, entry) in enumerate(stencil)
        selector, _ = _entry_selector_and_coeff(entry, i)
        axis_raw = get(selector, "axis", nothing)
        axis_raw isa AbstractString || throw(
            ArgumentError(
                "lower_stencil_to_replacement: arakawa entry $i 'selector.axis' " *
                    "must be a string",
            ),
        )
        axis = String(axis_raw)
        startswith(axis, "\$") || throw(
            ArgumentError(
                "lower_stencil_to_replacement: arakawa entry $i 'selector.axis' " *
                    "must be a '\$'-prefixed pattern variable (got '$axis')",
            ),
        )
        push!(seen, axis)
    end
    return sort!(collect(seen))
end

function _lower_arakawa_entry(
        selector::AbstractDict,
        coeff,
        operand_var::String,
        axis_vars::Vector{String},
        idx::Int,
    )
    sel_axis_raw = get(selector, "axis", nothing)
    sel_axis_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.axis' " *
                "must be a string",
        ),
    )
    sel_axis = String(sel_axis_raw)
    startswith(sel_axis, "\$") || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.axis' " *
                "must be a '\$'-prefixed pattern variable (got '$sel_axis')",
        ),
    )

    stagger_raw = get(selector, "stagger", nothing)
    stagger_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.stagger' " *
                "must be a string",
        ),
    )
    stagger = String(stagger_raw)
    stagger in _ARAKAWA_STAGGERS || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.stagger' " *
                "= '$stagger' is not one of $(_ARAKAWA_STAGGERS)",
        ),
    )

    offset_raw = get(selector, "offset", nothing)
    (offset_raw isa Integer && !(offset_raw isa Bool)) || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.offset' " *
                "must be an integer",
        ),
    )
    offset = Int(offset_raw)

    # Build the per-axis index slot list. The entry's `axis` slot carries
    # `$axis + offset` (or just `$axis` if offset == 0); all other slots
    # carry the bare axis pattern variable. Stagger info is intentionally
    # NOT in the AST — the binding-side runner recovers it from the rule
    # dict's preserved `stencil` field per SELECTOR_KINDS.md decision #16.
    index_args = Any[operand_var]
    sel_axis_in_axis_vars = false
    for ax in axis_vars
        if ax == sel_axis
            sel_axis_in_axis_vars = true
            push!(
                index_args,
                offset == 0 ? ax :
                    Dict{String, Any}("op" => "+", "args" => Any[ax, offset]),
            )
        else
            push!(index_args, ax)
        end
    end
    sel_axis_in_axis_vars || throw(
        ArgumentError(
            "lower_stencil_to_replacement: arakawa entry $idx 'selector.axis' " *
                "= '$sel_axis' is not in the stencil's collected axis set " *
                "$(axis_vars) (this should not happen)",
        ),
    )

    index_node = Dict{String, Any}("op" => "index", "args" => index_args)
    return Dict{String, Any}("op" => "*", "args" => Any[coeff, index_node])
end
