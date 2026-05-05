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

Selector kinds currently supported: `cartesian`. Other families
(`vertical`, `latlon`, `cubed_sphere`, `arakawa`, `indirect`,
`reduction`) raise `ArgumentError` until their lowering rows are added;
each new kind composes here as a separate dispatch branch.
"""

"""
    lower_stencil_to_replacement(rule::AbstractDict) -> Dict{String,Any}

Return a copy of `rule` with a `replacement` field inserted. Idempotent:
if `rule` already carries `replacement`, the input is returned via
shallow copy (no recomputation).

Stencil-form lowering produces

    combine_op(coeff_i * index(operand, axis_pat + offset_i)
               for entry i in stencil)

where:

- `combine_op` is the rule's top-level `combine` field (defaults to `"+"`).
- `operand` is the first `\$`-prefixed positional argument of `applies_to.args`.
- `axis_pat` is `applies_to.dim` (must itself be a `\$`-prefixed pattern variable).
- Each entry's `selector.axis` must equal `axis_pat`.

A stencil with a single entry lowers to the bare term — no `combine_op`
wrapper is emitted. An offset of `0` lowers to `index(operand, axis_pat)`
(the trivial `+ 0` is omitted).

Errors:

- `ArgumentError` if `rule` has neither `stencil` nor `replacement`.
- `ArgumentError` if `applies_to`, `applies_to.args`, or `applies_to.dim`
  is missing or malformed.
- `ArgumentError` if any selector kind is not yet supported.
- `ArgumentError` if any selector axis disagrees with `applies_to.dim`.
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
    axis_var = _stencil_axis_pattern_var(applies_to)

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

    terms = Any[]
    for (i, entry) in enumerate(stencil)
        push!(terms, _lower_stencil_entry(entry, operand_var, axis_var, i))
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

function _lower_stencil_entry(
        entry, operand_var::String, axis_var::String, idx::Int
    )
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

    kind = String(get(selector, "kind", ""))
    if kind == "cartesian"
        return _lower_cartesian_entry(selector, coeff, operand_var, axis_var, idx)
    end

    throw(
        ArgumentError(
            "lower_stencil_to_replacement: stencil entry $idx selector kind " *
                "'$kind' is not yet supported (currently supported: 'cartesian')",
        ),
    )
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
