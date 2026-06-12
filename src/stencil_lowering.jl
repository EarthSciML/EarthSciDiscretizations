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

- [`lower_stencil_to_replacement`](@ref) — single-rule lowerer to the
  scalar `replacement` AST form (legacy path; per-cell scalarization).
- [`lower_stencil_to_scheme`](@ref) — single-rule lowerer to the ESS §7
  scheme + `use:` rule document parts consumed natively by ESS
  `discretize` (scheme expansion → ArrayOp lift → tree-walk eval). This
  is the canonical-pipeline path; currently cartesian-only, mirroring
  the ESS scheme-expansion foundation (esm-j1u).
- [`lower_stencil_to_canonical_replacement`](@ref) — replacement
  expression with axis pattern variables substituted by canonical
  `\$target` component names (`i`, `j`, …; RFC §7.1). The form multi-axis
  (dim-less pattern) rules need to ride ESS `discretize` as plain
  `pattern` + `replacement` document rules.

Selector kinds currently supported by `lower_stencil_to_replacement`:
`cartesian`, `arakawa`, `latlon`, `cubed_sphere`, `vertical`. Other
families (`indirect`, `reduction`) raise `ArgumentError` until their
lowering rows are added; each new kind composes here as a separate
dispatch branch.
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

- **latlon** (`kind="latlon"`): multi-axis with literal geographic axis
  names (e.g. `"lon"`, `"lat"`). Unlike arakawa, the axis field is a
  literal string (not a `\$`-prefixed pattern variable). The sorted set
  of unique axis names defines the positional dimension order (`"lat"`
  sorts before `"lon"`). Each entry
  contributes `index(operand, lat_arg, lon_arg)` where the entry's axis
  slot carries `axis + offset` and other slots carry the bare literal axis
  name. No `stagger` field is required. `applies_to.dim` is present but
  not used — axis names are intrinsic to the stencil entries.

- **cubed_sphere** (`kind="cubed_sphere"`): multi-axis using a plural
  `selectors` array (not the singular `selector` key). Each stencil entry
  carries one selector per axis (e.g. `xi` and `eta`), each specifying an
  independent integer offset. The sorted set of all unique axis names
  defines the positional dimension order (`"eta"` sorts before `"xi"`).
  Each entry contributes
  `index(operand, eta_arg, xi_arg)` built from its own selector offsets.
  `applies_to.dim` is not required.

- **vertical** (`kind="vertical"`): single-axis with a required `stagger`
  field (`face_bottom`, `face_top`, `cell_center`). Uses `applies_to.dim`
  as the axis pattern variable (like cartesian). The stagger is NOT
  encoded in the lowered AST — the binding-side runner recovers it from
  the preserved `stencil` field to choose the correct physical array
  (analogous to arakawa stagger per decision #16).

A stencil with a single entry lowers to the bare term — no `combine_op`
wrapper is emitted. An offset of `0` lowers to `axis_pat` (the trivial
`+ 0` is omitted) at that axis slot.

Errors:

- `ArgumentError` if `rule` has neither `stencil` nor `replacement`.
- `ArgumentError` if `applies_to` or `applies_to.args` is missing or
  malformed.
- `ArgumentError` if entries declare conflicting selector kinds.
- `ArgumentError` if any selector kind is not yet supported.
- `ArgumentError` if a cartesian or vertical rule is missing
  `applies_to.dim` or has an entry whose axis disagrees with it.
- `ArgumentError` if an arakawa entry's `axis` is not a `\$`-prefixed
  pattern variable, or its `stagger` is not one of the four allowed
  symbols.
- `ArgumentError` if a latlon entry's `axis` is a `\$`-prefixed variable
  (must be a literal geographic axis name).
- `ArgumentError` if a cubed_sphere entry is missing its `selectors` array
  or a selector within it is missing its axis or offset.
- `ArgumentError` if a vertical entry's `stagger` is not one of
  `"face_bottom"`, `"face_top"`, `"cell_center"`.
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
    elseif family == "latlon"
        axis_names = _latlon_collect_axis_names(stencil)
        Any[_lower_latlon_entry_unpack(entry, operand_var, axis_names, i)
            for (i, entry) in enumerate(stencil)]
    elseif family == "cubed_sphere"
        axis_names = _cubed_sphere_collect_axis_names(stencil)
        Any[_lower_cubed_sphere_entry_unpack(entry, operand_var, axis_names, i)
            for (i, entry) in enumerate(stencil)]
    elseif family == "vertical"
        axis_var = _stencil_axis_pattern_var(applies_to)
        Any[_lower_vertical_entry_unpack(entry, operand_var, axis_var, i)
            for (i, entry) in enumerate(stencil)]
    else
        throw(
            ArgumentError(
                "lower_stencil_to_replacement: stencil selector kind " *
                    "'$family' is not yet supported (currently supported: " *
                    "'cartesian', 'arakawa', 'latlon', 'cubed_sphere', 'vertical')",
            ),
        )
    end

    replacement = length(terms) == 1 ? terms[1] :
                  Dict{String, Any}("op" => combine_op, "args" => terms)

    out["replacement"] = replacement
    return out
end

"""
    lower_stencil_to_scheme(name::AbstractString, rule::AbstractDict;
                            output::Union{Nothing,AbstractString}=nothing)
        -> (scheme::Dict{String,Any}, use_rule::Dict{String,Any})

Lower a stencil-form catalog rule to the two ESM document parts ESS's
canonical pipeline consumes natively (RFC §7.2.1, ESS scheme-expansion
foundation esm-j1u):

- `scheme` — the entry to place under the document's
  `discretizations.<name>` block. A whitelisted structural subset of the
  catalog spec: `applies_to`, `grid_family`, `stencil`, `combine`, and
  (when present) `accuracy` / `order`. Catalog-side extras such as a
  sibling `replacement` AST are dropped — ESS `parse_scheme` owns the
  schema from here.
- `use_rule` — the entry to append to the document's `rules` array:
  `{"name": name, "pattern": <applies_to>, "use": name}`. The pattern is
  the scheme's `applies_to` verbatim, so ESS `_check_applies_to` is
  satisfied by construction and the pattern variables bind one-to-one
  to the scheme operands.

This is a pure structural transform — no math is computed, no selector
is materialized; expansion to `coeff * index(...)` terms happens inside
ESS `expand_scheme`. Mirroring the ESS cartesian foundation, only
`grid_family == "cartesian"` rules whose selectors are all
`kind == "cartesian"` are lowerable; anything else raises
`ArgumentError` (the per-family extension lands when ESS grows the
corresponding selector dispatch — cubed-sphere panel esm-57f,
unstructured esm-bpr).

Multi-output rules (catalog extension ahead of the RFC: `stencil` is an
**object keyed by output name**, e.g. `ppm_reconstruction`'s
`q_left_edge` / `q_right_edge`) lower **one output at a time**: pass
`output=<name>` to select the entry list. Each output is an ordinary
single-axis cartesian stencil from ESS's point of view, so the emitted
scheme + `use:` rule drive the canonical pipeline unchanged; one ESM
document carries one output's scheme. Document-level multi-output
emission (a single rewrite producing every named output) awaits an RFC
§7 extension and is not expressible here.

Errors (`ArgumentError`): missing/empty `stencil`; malformed
`applies_to` (no `\$`-operand or non-`\$` `dim`); `grid_family` other
than `"cartesian"`; any selector with `kind != "cartesian"`, a
non-integer `offset`, or an `axis` disagreeing with `applies_to.dim`;
unsupported `combine`.
"""
function lower_stencil_to_scheme(name::AbstractString, rule::AbstractDict;
                                 output::Union{Nothing, AbstractString} = nothing)
    haskey(rule, "stencil") || throw(
        ArgumentError(
            "lower_stencil_to_scheme: rule has no 'stencil' array (replacement-form " *
                "rules go through the document 'rules' path directly)",
        ),
    )

    applies_to = get(rule, "applies_to", nothing)
    applies_to isa AbstractDict || throw(
        ArgumentError(
            "lower_stencil_to_scheme: rule missing 'applies_to' object",
        ),
    )
    _stencil_operand_pattern_var(applies_to)  # validates the $-operand exists
    axis_var = _stencil_axis_pattern_var(applies_to)

    grid_family = String(get(rule, "grid_family", ""))
    grid_family == "cartesian" || throw(
        ArgumentError(
            "lower_stencil_to_scheme: grid_family '$grid_family' is not lowerable " *
                "to an ESS scheme (ESS scheme expansion is cartesian-only; " *
                "cubed-sphere tracked at ESS/esm-57f, unstructured at ESS/esm-bpr)",
        ),
    )

    stencil_field = rule["stencil"]
    stencil = if stencil_field isa AbstractDict
        output === nothing && throw(
            ArgumentError(
                "lower_stencil_to_scheme: rule carries a multi-output stencil object " *
                    "(outputs: $(join(sort!(collect(String.(keys(stencil_field)))), ", "))); " *
                    "pass output=<name> to select one",
            ),
        )
        haskey(stencil_field, output) || throw(
            ArgumentError(
                "lower_stencil_to_scheme: no stencil output named '$output' " *
                    "(available: $(join(sort!(collect(String.(keys(stencil_field)))), ", ")))",
            ),
        )
        stencil_field[output]
    else
        output === nothing || throw(
            ArgumentError(
                "lower_stencil_to_scheme: output='$output' given but the rule's " *
                    "stencil is a flat (single-output) array",
            ),
        )
        stencil_field
    end
    stencil isa AbstractVector && !isempty(stencil) || throw(
        ArgumentError(
            "lower_stencil_to_scheme: 'stencil' must be a non-empty array",
        ),
    )

    combine_op = String(get(rule, "combine", "+"))
    combine_op in ("+", "*", "min", "max") || throw(
        ArgumentError(
            "lower_stencil_to_scheme: unsupported combine '$combine_op' " *
                "(expected one of '+', '*', 'min', 'max')",
        ),
    )

    lowered_stencil = Any[]
    for (i, entry) in enumerate(stencil)
        entry isa AbstractDict || throw(
            ArgumentError("lower_stencil_to_scheme: stencil entry $i is not an object"),
        )
        selector = get(entry, "selector", nothing)
        selector isa AbstractDict || throw(
            ArgumentError("lower_stencil_to_scheme: stencil entry $i has no 'selector' object"),
        )
        kind = String(get(selector, "kind", ""))
        kind == "cartesian" || throw(
            ArgumentError(
                "lower_stencil_to_scheme: stencil entry $i has selector kind '$kind' " *
                    "(ESS scheme expansion accepts 'cartesian' only)",
            ),
        )
        axis = get(selector, "axis", nothing)
        axis == axis_var || throw(
            ArgumentError(
                "lower_stencil_to_scheme: stencil entry $i selector axis '$axis' " *
                    "disagrees with applies_to.dim '$axis_var'",
            ),
        )
        offset = get(selector, "offset", nothing)
        offset isa Integer || throw(
            ArgumentError(
                "lower_stencil_to_scheme: stencil entry $i selector offset must be " *
                    "an integer (got $(repr(offset)))",
            ),
        )
        haskey(entry, "coeff") || throw(
            ArgumentError("lower_stencil_to_scheme: stencil entry $i has no 'coeff'"),
        )
        push!(lowered_stencil, Dict{String, Any}(
            "selector" => Dict{String, Any}(
                "kind"   => "cartesian",
                "axis"   => String(axis),
                "offset" => Int(offset),
            ),
            "coeff" => entry["coeff"],
        ))
    end

    scheme = Dict{String, Any}(
        "applies_to"  => applies_to,
        "grid_family" => "cartesian",
        "combine"     => combine_op,
        "stencil"     => lowered_stencil,
    )
    for passthrough in ("accuracy", "order")
        haskey(rule, passthrough) && (scheme[passthrough] = rule[passthrough])
    end

    use_rule = Dict{String, Any}(
        "name"    => String(name),
        "pattern" => applies_to,
        "use"     => String(name),
    )

    return scheme, use_rule
end

"""
    lower_stencil_to_canonical_replacement(rule::AbstractDict) -> Dict{String,Any}

Lower a stencil-form catalog rule to a replacement-expression AST whose
axis references use ESS's canonical `\$target` component names (`i`, `j`,
`k`, `l`, `m` — RFC §7.1 reserved keywords, position-based), instead of
the rule's axis pattern variables.

This is the form a multi-axis rule needs to ride ESS `discretize` as a
plain `pattern` + `replacement` document rule: a dim-less pattern like
`laplacian(\$u)` binds no axis variables, so the replacement cannot
reference `\$x`/`\$y` — but canonical component names are resolved
positionally by the arrayop lift regardless of the document's dimension
names (mirroring ESS's own 2D e2e rule authoring).

Pipeline: `lower_stencil_to_replacement` (idempotent) → strip the
`arrayop` wrapper if present → substitute each axis pattern variable
with its canonical component name. The substitution order is the
lowerer's documented positional dimension order (sorted unique axis
pattern variables — e.g. `\$x → i`, `\$y → j`). Returns the bare
replacement **expression** (no arrayop wrapper).

Restrictions (`ArgumentError`):

- Only collocated stencils: any `arakawa` selector with a `stagger`
  other than `"cell_center"` is rejected — canonical cell-component
  indexing into a face/corner-staggered array would silently read the
  wrong physical location.
- Only `\$`-pattern-variable axes (`cartesian`/`arakawa` kinds). Rules
  with literal axis names (`latlon`) or plural `selectors`
  (`cubed_sphere`) keep their existing lowering paths.
- At most 5 distinct axes (the canonical component alphabet).
"""
function lower_stencil_to_canonical_replacement(rule::AbstractDict)::Dict{String, Any}
    canonical_components = ("i", "j", "k", "l", "m")

    axis_vars = String[]
    stencil = get(rule, "stencil", nothing)
    if stencil isa AbstractVector
        seen = Set{String}()
        for (idx, entry) in enumerate(stencil)
            entry isa AbstractDict || continue
            selector = get(entry, "selector", nothing)
            selector isa AbstractDict || continue
            kind = String(get(selector, "kind", ""))
            if kind == "arakawa"
                stagger = String(get(selector, "stagger", ""))
                stagger == "cell_center" || throw(
                    ArgumentError(
                        "lower_stencil_to_canonical_replacement: stencil entry $idx has " *
                            "stagger '$stagger' — canonical-component form is only valid " *
                            "for collocated (cell_center) stencils",
                    ),
                )
            end
            axis = get(selector, "axis", nothing)
            if axis isa AbstractString && startswith(String(axis), "\$")
                a = String(axis)
                a in seen || (push!(axis_vars, a); push!(seen, a))
            elseif axis !== nothing
                throw(
                    ArgumentError(
                        "lower_stencil_to_canonical_replacement: stencil entry $idx axis " *
                            "$(repr(axis)) is not a '\$'-prefixed pattern variable",
                    ),
                )
            end
        end
    else
        # Replacement-form rule: the only axis variable is applies_to.dim.
        applies_to = get(rule, "applies_to", nothing)
        if applies_to isa AbstractDict
            dim = get(applies_to, "dim", nothing)
            if dim isa AbstractString && startswith(String(dim), "\$")
                push!(axis_vars, String(dim))
            end
        end
    end
    sort!(axis_vars)
    length(axis_vars) <= length(canonical_components) || throw(
        ArgumentError(
            "lower_stencil_to_canonical_replacement: $(length(axis_vars)) distinct axes " *
                "exceed the canonical component alphabet (i, j, k, l, m)",
        ),
    )
    subst = Dict{String, String}(
        v => canonical_components[k] for (k, v) in enumerate(axis_vars))

    lowered = lower_stencil_to_replacement(rule)
    repl = lowered["replacement"]
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl

    return _substitute_axis_vars(expr, subst)
end

# Replace axis pattern-variable string leaves with their canonical component
# names; every other node passes through structurally unchanged.
function _substitute_axis_vars(expr, subst::Dict{String, String})
    if expr isa AbstractString
        return get(subst, String(expr), String(expr))
    end
    expr isa AbstractDict || return expr
    out = Dict{String, Any}()
    for (k, v) in expr
        key = String(k)
        if key == "args" && v isa AbstractVector
            out[key] = Any[_substitute_axis_vars(a, subst) for a in v]
        else
            out[key] = v
        end
    end
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

Inspect every entry's selector kind and return the (single) family name
shared across the stencil. All entries must agree; mixed-kind stencils
are not supported.

Handles both the singular `selector` key (cartesian, arakawa, latlon,
vertical) and the plural `selectors` array (cubed_sphere).
"""
function _stencil_family(stencil::AbstractVector)::String
    family = ""
    for (i, entry) in enumerate(stencil)
        entry isa AbstractDict || throw(
            ArgumentError(
                "lower_stencil_to_replacement: stencil entry $i must be an object",
            ),
        )
        kind = if haskey(entry, "selectors")
            sels = entry["selectors"]
            (sels isa AbstractVector && !isempty(sels)) || throw(
                ArgumentError(
                    "lower_stencil_to_replacement: stencil entry $i 'selectors' " *
                        "must be a non-empty array",
                ),
            )
            first_sel = sels[1]
            first_sel isa AbstractDict || throw(
                ArgumentError(
                    "lower_stencil_to_replacement: stencil entry $i 'selectors[1]' " *
                        "must be an object",
                ),
            )
            String(get(first_sel, "kind", ""))
        else
            selector = get(entry, "selector", nothing)
            selector isa AbstractDict || throw(
                ArgumentError(
                    "lower_stencil_to_replacement: stencil entry $i missing " *
                        "'selector' object (or 'selectors' array for multi-axis families)",
                ),
            )
            String(get(selector, "kind", ""))
        end
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

# -----------------------------------------------------------------------------
# Latlon lowering
# -----------------------------------------------------------------------------

"""
    _latlon_collect_axis_names(stencil) -> Vector{String}

Walk every entry's `selector.axis` and return the sorted set of unique
literal geographic axis names (e.g. `"lon"`, `"lat"`). Unlike arakawa,
latlon axes are NOT `\$`-prefixed pattern variables. `\$`-prefixed values
are rejected.
"""
function _latlon_collect_axis_names(stencil::AbstractVector)::Vector{String}
    seen = Set{String}()
    for (i, entry) in enumerate(stencil)
        selector, _ = _entry_selector_and_coeff(entry, i)
        axis_raw = get(selector, "axis", nothing)
        axis_raw isa AbstractString || throw(
            ArgumentError(
                "lower_stencil_to_replacement: latlon entry $i 'selector.axis' " *
                    "must be a string",
            ),
        )
        axis = String(axis_raw)
        startswith(axis, "\$") && throw(
            ArgumentError(
                "lower_stencil_to_replacement: latlon entry $i 'selector.axis' " *
                    "must be a literal geographic axis name, not a '\$'-prefixed " *
                    "pattern variable (got '$axis')",
            ),
        )
        push!(seen, axis)
    end
    return sort!(collect(seen))
end

function _lower_latlon_entry_unpack(
        entry, operand_var::String, axis_names::Vector{String}, idx::Int)
    selector, coeff = _entry_selector_and_coeff(entry, idx)
    return _lower_latlon_entry(selector, coeff, operand_var, axis_names, idx)
end

function _lower_latlon_entry(
        selector::AbstractDict,
        coeff,
        operand_var::String,
        axis_names::Vector{String},
        idx::Int,
    )
    sel_axis_raw = get(selector, "axis", nothing)
    sel_axis_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: latlon entry $idx 'selector.axis' " *
                "must be a string",
        ),
    )
    sel_axis = String(sel_axis_raw)
    startswith(sel_axis, "\$") && throw(
        ArgumentError(
            "lower_stencil_to_replacement: latlon entry $idx 'selector.axis' " *
                "must be a literal geographic axis name, not a '\$'-prefixed " *
                "pattern variable (got '$sel_axis')",
        ),
    )

    offset_raw = get(selector, "offset", nothing)
    (offset_raw isa Integer && !(offset_raw isa Bool)) || throw(
        ArgumentError(
            "lower_stencil_to_replacement: latlon entry $idx 'selector.offset' " *
                "must be an integer",
        ),
    )
    offset = Int(offset_raw)

    # Same multi-axis pattern as arakawa but with literal axis names.
    # Entry's axis slot carries `axis + offset` (or bare `axis` if offset == 0);
    # all other axis slots carry the bare literal axis name.
    index_args = Any[operand_var]
    sel_axis_found = false
    for ax in axis_names
        if ax == sel_axis
            sel_axis_found = true
            push!(
                index_args,
                offset == 0 ? ax :
                    Dict{String, Any}("op" => "+", "args" => Any[ax, offset]),
            )
        else
            push!(index_args, ax)
        end
    end
    sel_axis_found || throw(
        ArgumentError(
            "lower_stencil_to_replacement: latlon entry $idx 'selector.axis' " *
                "= '$sel_axis' is not in the stencil's collected axis set " *
                "$(axis_names) (this should not happen)",
        ),
    )

    index_node = Dict{String, Any}("op" => "index", "args" => index_args)
    return Dict{String, Any}("op" => "*", "args" => Any[coeff, index_node])
end

# -----------------------------------------------------------------------------
# Cubed-sphere lowering
# -----------------------------------------------------------------------------

"""
    _cubed_sphere_entry_selectors_and_coeff(entry, idx) -> (AbstractVector, coeff)

For cubed_sphere entries that use the plural `selectors` key (an array of
per-axis selector objects), extract the selectors array and the coefficient.
"""
function _cubed_sphere_entry_selectors_and_coeff(entry, idx::Int)
    entry isa AbstractDict || throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx must be an object",
        ),
    )
    selectors = get(entry, "selectors", nothing)
    selectors isa AbstractVector || throw(
        ArgumentError(
            "lower_stencil_to_replacement: cubed_sphere stencil entry $idx " *
                "missing 'selectors' array",
        ),
    )
    coeff = get(entry, "coeff", nothing)
    coeff === nothing && throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx missing 'coeff'",
        ),
    )
    return selectors, coeff
end

"""
    _cubed_sphere_collect_axis_names(stencil) -> Vector{String}

Walk every entry's `selectors` array and return the sorted set of unique
literal axis names (e.g. `"eta"`, `"xi"`).
"""
function _cubed_sphere_collect_axis_names(stencil::AbstractVector)::Vector{String}
    seen = Set{String}()
    for (i, entry) in enumerate(stencil)
        selectors, _ = _cubed_sphere_entry_selectors_and_coeff(entry, i)
        for (j, sel) in enumerate(selectors)
            sel isa AbstractDict || throw(
                ArgumentError(
                    "lower_stencil_to_replacement: cubed_sphere entry $i " *
                        "selector $j must be an object",
                ),
            )
            axis_raw = get(sel, "axis", nothing)
            axis_raw isa AbstractString || throw(
                ArgumentError(
                    "lower_stencil_to_replacement: cubed_sphere entry $i " *
                        "selector $j 'axis' must be a string",
                ),
            )
            push!(seen, String(axis_raw))
        end
    end
    return sort!(collect(seen))
end

function _lower_cubed_sphere_entry_unpack(
        entry, operand_var::String, axis_names::Vector{String}, idx::Int)
    selectors, coeff = _cubed_sphere_entry_selectors_and_coeff(entry, idx)
    return _lower_cubed_sphere_entry(selectors, coeff, operand_var, axis_names, idx)
end

function _lower_cubed_sphere_entry(
        selectors::AbstractVector,
        coeff,
        operand_var::String,
        axis_names::Vector{String},
        idx::Int,
    )
    # Build a map from axis name -> offset for this entry from its selectors array.
    axis_offset = Dict{String, Int}()
    for (j, sel) in enumerate(selectors)
        sel isa AbstractDict || throw(
            ArgumentError(
                "lower_stencil_to_replacement: cubed_sphere entry $idx " *
                    "selector $j must be an object",
            ),
        )
        axis_raw = get(sel, "axis", nothing)
        axis_raw isa AbstractString || throw(
            ArgumentError(
                "lower_stencil_to_replacement: cubed_sphere entry $idx " *
                    "selector $j 'axis' must be a string",
            ),
        )
        axis = String(axis_raw)
        offset_raw = get(sel, "offset", nothing)
        (offset_raw isa Integer && !(offset_raw isa Bool)) || throw(
            ArgumentError(
                "lower_stencil_to_replacement: cubed_sphere entry $idx " *
                    "selector $j 'offset' must be an integer",
            ),
        )
        axis_offset[axis] = Int(offset_raw)
    end

    # Build index args in sorted axis order; every axis must appear in this entry.
    index_args = Any[operand_var]
    for ax in axis_names
        offset = get(axis_offset, ax, nothing)
        offset === nothing && throw(
            ArgumentError(
                "lower_stencil_to_replacement: cubed_sphere entry $idx " *
                    "missing selector for axis '$ax' (expected in every entry)",
            ),
        )
        push!(
            index_args,
            offset == 0 ? ax :
                Dict{String, Any}("op" => "+", "args" => Any[ax, offset]),
        )
    end

    index_node = Dict{String, Any}("op" => "index", "args" => index_args)
    return Dict{String, Any}("op" => "*", "args" => Any[coeff, index_node])
end

# -----------------------------------------------------------------------------
# Vertical lowering
# -----------------------------------------------------------------------------

const _VERTICAL_STAGGERS = ("face_bottom", "face_top", "cell_center")

function _lower_vertical_entry_unpack(
        entry, operand_var::String, axis_var::String, idx::Int)
    selector, coeff = _entry_selector_and_coeff(entry, idx)
    return _lower_vertical_entry(selector, coeff, operand_var, axis_var, idx)
end

function _lower_vertical_entry(
        selector::AbstractDict,
        coeff,
        operand_var::String,
        axis_var::String,
        idx::Int,
    )
    sel_axis_raw = get(selector, "axis", nothing)
    sel_axis_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: vertical entry $idx 'selector.axis' " *
                "must be a string",
        ),
    )
    sel_axis = String(sel_axis_raw)
    if sel_axis != axis_var
        throw(
            ArgumentError(
                "lower_stencil_to_replacement: vertical entry $idx 'selector.axis' " *
                    "= '$sel_axis' disagrees with 'applies_to.dim' = '$axis_var'",
            ),
        )
    end

    stagger_raw = get(selector, "stagger", nothing)
    stagger_raw isa AbstractString || throw(
        ArgumentError(
            "lower_stencil_to_replacement: vertical entry $idx 'selector.stagger' " *
                "must be a string",
        ),
    )
    stagger = String(stagger_raw)
    stagger in _VERTICAL_STAGGERS || throw(
        ArgumentError(
            "lower_stencil_to_replacement: vertical entry $idx 'selector.stagger' " *
                "= '$stagger' is not one of $(_VERTICAL_STAGGERS)",
        ),
    )

    offset_raw = get(selector, "offset", nothing)
    (offset_raw isa Integer && !(offset_raw isa Bool)) || throw(
        ArgumentError(
            "lower_stencil_to_replacement: vertical entry $idx 'selector.offset' " *
                "must be an integer",
        ),
    )
    offset = Int(offset_raw)

    # Stagger is NOT encoded in the lowered AST — binding-side recovers it from
    # the preserved `stencil` field (same convention as arakawa, decision #16).
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
