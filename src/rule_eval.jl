"""
    rule_eval.jl — thin passthrough to the EarthSciSerialization AST evaluator.

ESD rule authors express stencil coefficients, edge interpolants, and reconstructions
directly as ExpressionNode ASTs inside rule JSON files (ESS §7 / §9.2). This module
exposes a single Julia entry point that delegates to the ESS evaluator instead of
carrying a shadow evaluator inside ESD.

Public surface:

- `eval_coeff(node, bindings)` — evaluate a JSON-decoded ExpressionNode (Number,
  String, or `Dict{String,Any}` with `"op"`/`"args"`/optional `"wrt"`/`"dim"`)
  against a `Dict{String,Float64}` of variable bindings. Returns `Float64`.

The heavy lifting — operator semantics, domain errors, unbound-variable detection —
lives in `EarthSciSerialization`. This file exists so the rest of ESD imports ONE
thing (`eval_coeff`) regardless of how ESS's public surface evolves.
"""

import EarthSciSerialization
const _ESS = EarthSciSerialization

"""
    eval_coeff(node, bindings::Dict{String,Float64}) -> Float64

Evaluate a JSON-decoded rule-AST coefficient node using the ESS tree-walk evaluator.

`node` is what `JSON.parse` produces from a rule file's `"coeff"`, `"stencil"` entry,
or any nested `{"op": ..., "args": ...}` block: a `Number`, a `String` (variable
name), or an `AbstractDict` with at minimum `"op"` and `"args"` keys.

`bindings` maps variable names to concrete floats (e.g. `Dict("dx" => 1/64)`).

Throws `EarthSciSerialization.UnboundVariableError` if a variable appears in the
AST but not in `bindings`, and propagates any numerical errors (`DomainError`,
`DivideError`) raised by the ESS evaluator.
"""
function eval_coeff(node, bindings::Dict{String,Float64})::Float64
    expr = _ESS.parse_expression(node)
    return _ESS.evaluate(expr, bindings)
end
