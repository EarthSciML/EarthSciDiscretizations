# Operators

## The closed authoring vocabulary

A discretization rule's `replacement` MUST be expressible in the §4.2 op
vocabulary alone. That vocabulary is *closed*: every binding (Julia,
Python, future host) walks the same AST, and no rule introduces new
ops. The vocabulary covers everything an FV scheme needs — including
limiters, reconstructions, and weighted-stencil schemes — without any
scheme-specific dispatch.

The full list, summarized from `esm-spec.md` §4.2 (consult that section
for the authoritative form):

### PDE operators (`applies_to.op` and tendency expressions)

| Op | Fields | Meaning |
|---|---|---|
| `D` | `wrt` | Time derivative ∂/∂t |
| `grad` | `dim` | Spatial gradient ∂/∂x |
| `div` | | Divergence ∇· |
| `laplacian` | | Laplacian ∇² |

These four (plus pointwise math, below) are the PDE-operator alphabet
that authors pattern-match against. A rule's `applies_to.op` MUST be
one of them.

### Pointwise math

Arithmetic: `+` (n-ary), `-` (unary or binary), `*` (n-ary), `/`
(binary), `^` (binary).

Elementary functions: `exp`, `log`, `log10`, `sqrt`, `abs`, `sign`,
`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `min` (n-ary,
≥2 args), `max` (n-ary, ≥2 args), `floor`, `ceil`.

`min` / `max` are the canonical encoding of clamp / clip / limiter
primitives — `clamp(x, lo, hi)` is `min(hi, max(lo, x))`. Reviewers
reject `fn` nodes that re-implement these in disguise.

### Conditionals

`ifelse [cond, then, else]`, `>`, `<`, `>=`, `<=`, `==`, `!=`, `and`,
`or`, `not`. Conditional logic in lowerings (e.g. flux-limiter switches)
goes through `ifelse`; there is no `if`-statement primitive.

### Inline constants and registry calls

`const` (literal value), `enum` (file-local symbolic name), `fn`
(invocation of a closed-registry function — `name` is a dotted path).
Most FV rules need none of these; mention them only when shipping a
table of coefficients inline (`const`) or a registry-level helper.

### Array / tensor ops (`replacement` body)

| Op | Required fields | Meaning |
|---|---|---|
| `arrayop` | `output_idx`, `expr` | Generalized Einstein-notation tensor expression with implicit reductions over non-output indices. (§4.3.1.) |
| `makearray` | `regions`, `values` | Block assembly from overlapping sub-region assignments; later regions overwrite earlier ones. (§4.3.2.) |
| `index` | — | Element / sub-array access. `args[0]` is the array; `args[1..]` are index expressions. (§4.3.3.) |
| `broadcast` | `fn` | Element-wise application of scalar op `fn`. (§4.3.4.) |
| `reshape` | `shape` | Reshape the operand to a target shape. (§4.3.5.) |
| `transpose` | optional `perm` | Axis permutation. (§4.3.5.) |
| `concat` | `axis` | Concatenate operand arrays along `axis`. (§4.3.5.) |

`arrayop` is the workhorse: every closed-AST FV lowering in the
catalog ultimately resolves to one or more `arrayop` nodes whose body
combines `index`, arithmetic, and (for nonlinear schemes) `ifelse` and
`broadcast`.

## What the vocabulary excludes

Rules MUST NOT introduce off-spec match keys or scheme-named ops in
either `applies_to.op` or anywhere inside `replacement`. The following
names are forbidden:

| Forbidden as `op` | Right answer |
|---|---|
| `advect` | match `D(q, wrt=t)` (or whichever PDE op the scheme discretizes) and lower the advective tendency as an `arrayop`. |
| `reconstruct` | the reconstruction is the body of the `arrayop` — express the polynomial in `+`, `*`, `index`. |
| `flux` | flux assembly is an `arrayop` whose `output_idx` ranges over edges and whose body is the numerical-flux formula in pointwise math + `ifelse`. |
| `limit` / `limiter` | encode the limiter as `min` / `max` / `ifelse` directly. |
| `bc:*` | boundary handling is declared on the domain (`esm-spec.md` §11.5) and applied by downstream BC rules; no `bc:*` op exists in the lowered AST. |

The rule of thumb: if the only justification for a new op name is "this
scheme calls it that", it is the wrong move. Express the math in the
existing alphabet. If a scheme genuinely cannot be expressed — bring it
to the spec authors before authoring the rule.

The companion prohibition is on **scheme-specific kernels in any host
language**. Rules ship JSON, not Julia or Python. The ESS rule engine's
AST walker implements the §4.2 vocabulary once; rules borrow that
implementation by composition. ESD carries no rule evaluator of its own
(`src/rule_eval.jl::eval_coeff` is a thin passthrough to the ESS
tree-walk evaluator).

## Worked sanity check

The closed-AST lowerings produced by the catalog satisfy the same basic
invariants any operator on arrays must. The catalog rules are evaluated
through the ESS rule engine via `build_ode_problem` (see
[Getting started: solve a PDE](@ref)) — there are
no named Julia operator callables to invoke directly. The example below
builds a real catalog rule (`centered_2nd_deriv_uniform`, the discrete
Laplacian) into an `ODEProblem` and confirms the basic invariant that
the discrete operator annihilates a constant field.

### A constant field has zero discrete second derivative

```@example ops
using EarthSciDiscretizations

repo = dirname(dirname(pathof(EarthSciDiscretizations)))
esm  = joinpath(repo, "test", "fixtures", "diffusion_1d.esm")
gdd  = joinpath(repo, "discretizations", "gdd", "cartesian_1d_periodic_n16.gdd.json")

prob, var_map = build_ode_problem(esm; grid_ref = gdd)

# Constant field over all 16 cells.
u  = fill(5.0, length(prob.u0))
du = similar(u)
prob.f(du, u, prob.p, 0.0)

println("Discrete ∂²u/∂x² of a constant field:")
println("  Size: $(size(du))")
println("  Max absolute value: $(maximum(abs.(du)))")
```

## Where to read more

- [Finite-Volume Method](@ref) — how a rule's pattern match and closed
  `arrayop` replacement encode an FV operator.
- [Tutorial: Authoring a rule](@ref) — end-to-end walkthrough.
- `esm-spec.md` §4.2 (operator vocabulary), §4.3 (array semantics).
