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
per-rule-shape logic.

Public surface (esd-t4h):

- [`lower_stencil_to_scheme`](@ref) — single-rule lowerer to the ESS §7
  scheme + `use:` rule document parts consumed natively by ESS
  `discretize` (scheme expansion → ArrayOp lift → tree-walk eval). This
  is the canonical-pipeline path; currently cartesian-only, mirroring
  the ESS scheme-expansion foundation (esm-j1u).

`lower_stencil_to_replacement` was retired in esd-t4h: all latlon/vertical
catalog rules now carry authored `replacement` ASTs, cubed_sphere stencil
rules are spec-only (EINSUM-8 territory), and arakawa rules were retired in
EINSUM-4 (esd-770). Only `lower_stencil_to_scheme` remains for the
cartesian/unstructured path.
"""

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
ESS `expand_scheme`. Supported `grid_family` values:

- **`"cartesian"`** — selectors must all be `kind == "cartesian"`;
  validates `axis` against `applies_to.dim` and `offset` as integer.
- **`"unstructured"`** — selectors must all be `kind == "reduction"` or
  `kind == "indirect"`; selector fields are passed through verbatim
  (`table`, `count_expr`, `k_bound`, `combine`, `index_expr`) for ESS
  `parse_scheme` to validate (ess-t0z). `applies_to.dim` need not be a
  `\$`-prefixed pattern variable.

Other `grid_family` values raise `ArgumentError`; latlon/vertical rules carry
authored `replacement` ASTs (esd-t4h) and flow through the document `rules`
path directly.

Multi-output rules (catalog extension: `stencil` is an **object keyed
by output name**, e.g. `ppm_reconstruction`'s `q_left_edge` /
`q_right_edge`) lower **one output at a time**: pass `output=<name>` to
select the entry list. Each output is an ordinary single-axis cartesian
stencil from ESS's point of view, so the emitted scheme + `use:` rule
drive the canonical pipeline unchanged; one ESM document carries one
output's scheme. This is the path the Layer-B convergence runner uses.

Document-level multi-output emission — a single ESM document carrying a
`kind: "multi_output_stencil"` scheme that ESS expands to observed
equations for every named output — is now expressible via the ESS
`multi_output_stencil` extension (ess-7hb/ess-ebe). The Layer-A
canonical fixture for `ppm_reconstruction` demonstrates this path: its
`input.esm` carries `kind: "multi_output_stencil"` with `primary:
"q_right_edge"`, and `discretize` expands it to observed arrayop
equations for both `q_left_edge__q__x` and `q_right_edge__q__x`. This
function (`lower_stencil_to_scheme`) continues to emit single-output
schemes for Layer-B since the convergence runner iterates each output
independently.

Errors (`ArgumentError`): missing/empty `stencil`; malformed
`applies_to` (no `\$`-operand); `grid_family` not in `"cartesian"` or
`"unstructured"`; for cartesian: non-`\$` `dim`, selector `kind !=
"cartesian"`, non-integer `offset`, or `axis` disagreeing with
`applies_to.dim`; for unstructured: selector `kind` not in
`"reduction"` / `"indirect"`; unsupported `combine`.
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

    grid_family = String(get(rule, "grid_family", ""))
    grid_family in ("cartesian", "unstructured") || throw(
        ArgumentError(
            "lower_stencil_to_scheme: grid_family '$grid_family' is not lowerable " *
                "to an ESS scheme (supported: 'cartesian', 'unstructured'; latlon/vertical " *
                "rules carry authored replacement ASTs since esd-t4h)",
        ),
    )

    axis_var = grid_family == "cartesian" ? _stencil_axis_pattern_var(applies_to) : ""

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
        haskey(entry, "coeff") || throw(
            ArgumentError("lower_stencil_to_scheme: stencil entry $i has no 'coeff'"),
        )
        kind = String(get(selector, "kind", ""))
        if grid_family == "cartesian"
            kind == "cartesian" || throw(
                ArgumentError(
                    "lower_stencil_to_scheme: stencil entry $i has selector kind '$kind' " *
                        "(ESS scheme expansion accepts 'cartesian' only for cartesian grid family)",
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
            push!(lowered_stencil, Dict{String, Any}(
                "selector" => Dict{String, Any}(
                    "kind"   => "cartesian",
                    "axis"   => String(axis),
                    "offset" => Int(offset),
                ),
                "coeff" => entry["coeff"],
            ))
        else  # grid_family == "unstructured"
            kind in ("reduction", "indirect") || throw(
                ArgumentError(
                    "lower_stencil_to_scheme: stencil entry $i has selector kind '$kind' " *
                        "(unstructured grid_family accepts 'reduction' and 'indirect' only)",
                ),
            )
            push!(lowered_stencil, Dict{String, Any}(
                "selector" => Dict{String, Any}(String(k) => v for (k, v) in selector),
                "coeff"    => entry["coeff"],
            ))
        end
    end

    scheme = Dict{String, Any}(
        "applies_to"  => applies_to,
        "grid_family" => grid_family,
        "combine"     => combine_op,
        "stencil"     => lowered_stencil,
    )
    for passthrough in ("accuracy", "order")
        haskey(rule, passthrough) && (scheme[passthrough] = rule[passthrough])
    end
    if grid_family == "unstructured"
        for passthrough in ("requires_locations", "emits_location")
            haskey(rule, passthrough) && (scheme[passthrough] = rule[passthrough])
        end
    end

    use_rule = Dict{String, Any}(
        "name"    => String(name),
        "pattern" => applies_to,
        "use"     => String(name),
    )

    return scheme, use_rule
end

# -----------------------------------------------------------------------------
# Internal helpers shared by lower_stencil_to_scheme
# -----------------------------------------------------------------------------

function _stencil_operand_pattern_var(applies_to::AbstractDict)::String
    args = get(applies_to, "args", nothing)
    args isa AbstractVector || throw(
        ArgumentError(
            "stencil lowering: 'applies_to.args' must be an array",
        ),
    )
    for a in args
        if a isa AbstractString && startswith(String(a), "\$")
            return String(a)
        end
    end
    throw(
        ArgumentError(
            "stencil lowering: 'applies_to.args' must contain at " *
                "least one '\$'-prefixed pattern variable (operand)",
        ),
    )
end

function _stencil_axis_pattern_var(applies_to::AbstractDict)::String
    dim = get(applies_to, "dim", nothing)
    dim isa AbstractString || throw(
        ArgumentError(
            "stencil lowering: 'applies_to.dim' must be a string",
        ),
    )
    s = String(dim)
    startswith(s, "\$") || throw(
        ArgumentError(
            "stencil lowering: 'applies_to.dim' must be a " *
                "'\$'-prefixed pattern variable (got '$s')",
        ),
    )
    return s
end
