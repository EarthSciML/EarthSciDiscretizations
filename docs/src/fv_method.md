# Finite-Volume Method

## Overview

EarthSciDiscretizations.jl supplies finite-volume (FV) discretization rules
that the EarthSciSerialization (ESS) rewriter applies to PDE equation trees.
A rule is a small JSON file with two halves:

1. **`applies_to`** — a *pattern match* against a §4.2 PDE operator
   (`D`, `grad`, `div`, `laplacian`, plus the pointwise-math vocabulary).
2. **`replacement`** — a **closed `arrayop` lowering** in the §4.2 op
   vocabulary only (`arrayop`, `broadcast`, `index`, `ifelse`, `+`, `-`,
   `*`, `/`, `^`, `sqrt`, `min`, `max`, …).

Authors do **not** invent scheme-specific match keys (no `advect`,
`reconstruct`, `flux`, `limit` — those names are forbidden as `applies_to.op`
values) and do **not** ship scheme-specific kernels in any host language.
Every binding (Julia, Python, future host) executes the same closed AST by
walking the §4.2 ops it already understands.

The integral conservation law that motivates the FV method survives this
move — it now lives in the structure of the `arrayop` body rather than in a
named coefficient table.

## The integral form

For a vector field $\mathbf{F}$ over a cell with area $A$ and boundary
$\partial A$, the divergence theorem reads

```math
\int_A \nabla \cdot \mathbf{F} \, dA = \oint_{\partial A} \mathbf{F} \cdot \hat{n} \, ds.
```

A rule that lowers `div(F)` is therefore expected to produce an `arrayop`
whose body sums signed edge fluxes scaled by edge lengths and the inverse
cell area. The rule does not need a custom dispatcher to express this: the
`+` and `*` and `index` ops in §4.2 are sufficient.

## Worked example: `centered_2nd_uniform`

The smallest rule in the catalog —
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}),
landed in commit `7b26ffd` — is the canonical exemplar of the closed-AST
lowering pattern. It targets `grad(u, dim=x)` on a uniform Cartesian axis.

### (a) The PDE operator

The continuous operator is the partial derivative

```math
\left(\frac{\partial u}{\partial x}\right)(x).
```

In §4.2 it is encoded as `{"op": "grad", "args": ["u"], "dim": "x"}`.

### (b) The pattern match

The rule's `applies_to` clause matches that op verbatim, with metavariables
`$u` and `$x` for the operand and axis:

```json
"applies_to": {
  "op":   "grad",
  "args": ["$u"],
  "dim":  "$x"
}
```

The matcher fires anywhere in the equation tree where a `grad` node has the
declared shape; `$u` and `$x` bind to the concrete field name and axis at
rewrite time.

### (c) The closed `arrayop` replacement

The replacement is a single `arrayop` node whose body is a centered
two-point difference:

```json
"replacement": {
  "op": "arrayop",
  "output_idx": ["$x"],
  "expr": {
    "op": "/",
    "args": [
      { "op": "-", "args": [
        { "op": "index", "args": ["$u", { "op": "+", "args": ["$x", 1] }] },
        { "op": "index", "args": ["$u", { "op": "-", "args": ["$x", 1] }] }
      ]},
      { "op": "*", "args": [2, "dx"] }
    ]
  },
  "args": ["$u"]
}
```

After substitution (`$u → u`, `$x → i`), the lowering at each interior
cell reduces to

```math
\left(\frac{\partial u}{\partial x}\right)_i \approx \frac{u_{i+1} - u_{i-1}}{2\,\Delta x}.
```

There is no `coeff` table, no `stencil[]` array, no per-host kernel: the
arithmetic that yields the second-order centered difference is visible
directly in the AST.

## Boundary conditions live on the domain

A rule's replacement is the **interior** closed form. Boundary handling is
declared once per field on the domain's `boundary_conditions` block
(`esm-spec.md` §11.5: `periodic`, `dirichlet` / `constant`, `neumann` /
`zero_gradient`, `robin`) and applied as **downstream rewrite rules over
concrete indices**. The lowered AST does not contain `bc:*` ops.

| Domain BC | Index transformation applied to `$u[$x ± 1]` |
|---|---|
| `periodic` | wrap-around: `mod($x ± 1 + N, N)` (see [`periodic_bc`]({{< ref "/rules/periodic_bc" >}})) |
| `dirichlet` / `constant` | boundary cell reads the prescribed value |
| `neumann` / `zero_gradient` | mirror in-range neighbor (clamp the index) |
| `robin` | mixed coefficient row at the boundary |

This is the same separation the centered exemplar uses: the rule itself
stays BC-agnostic, and a downstream BC rule fires at the boundary cells.
Authors of new rules should not embed BC logic in their lowering.

## What "discrete FV operator" means in this framing

When this guide talks about a discrete divergence, gradient, or Laplacian,
it is shorthand for "the closed `arrayop` lowering produced when the
matching rule fires against the corresponding §4.2 op". For example:

- `grad(u, dim=x)` on a uniform Cartesian axis → `centered_2nd_uniform`
  → centered two-point difference (above).
- `div(F^ξ, F^η)` on a cubed-sphere C-grid → a rule that lowers to an
  `arrayop` summing edge-flux-times-edge-length and dividing by cell area.
- `laplacian(φ)` on the cubed sphere → a rule whose `arrayop` encodes the
  5-point orthogonal stencil plus the cross-metric correction terms (the
  geometry lives in the body of the AST, not in a separate dispatch table).

In every case the *math* — the integral form, the metric tensor, the
truncation error — is the motivation for the AST shape; the *authoring
artifact* is just the AST.

## Where to read more

- [Operators](@ref) — the §4.2 op vocabulary you may use inside a `replacement`.
- [Tutorial: Authoring a rule](@ref) — end-to-end walkthrough.
- `esm-spec.md` §4.2 (operator vocabulary), §4.3 (array semantics),
  §11.5 (BC types) for the definitive specification.

## C-grid staggering

The Arakawa C-grid staggering places different variables at different
locations within each cell. Index symbols in an `arrayop` lowering are
local to the rule, but the runtime resolves their lengths from the
operand shapes — which in turn come from the staggering:

| Location | Symbol | Grid Size | Description |
|:---------|:-------|:----------|:------------|
| `CellCenter` | $(i, j)$ | $(N_c, N_c)$ | Scalar fields (tracer, pressure, temperature) |
| `UEdge` | $(i+1/2, j)$ | $(N_c+1, N_c)$ | Normal velocity component in $\xi$-direction |
| `VEdge` | $(i, j+1/2)$ | $(N_c, N_c+1)$ | Normal velocity component in $\eta$-direction |
| `Corner` | $(i+1/2, j+1/2)$ | $(N_c+1, N_c+1)$ | Vorticity, stream function |

A rule that produces an edge-quantity output gives an `output_idx` whose
range matches the corresponding edge-staggered grid; one that consumes an
edge quantity uses `index` into an array shaped that way.

## Ghost cells

Inter-panel communication is handled through ghost cells. Each panel is
padded with $N_g$ ghost layers on each side, filled from neighboring
panels using the connectivity table and index transformations. This
happens at the *grid* level, before any rule fires; rule lowerings see
the ghosted arrays as ordinary inputs.

## References

The finite-volume methods this package targets are based on the following
foundational works. Their algorithms motivate the *shape* of the AST
lowerings; nothing in the rule files reproduces a host-language
implementation of them.

- Lin, S.-J. and R. B. Rood (1996). "Multidimensional Flux-Form Semi-Lagrangian Transport Schemes." *Monthly Weather Review*, 124(9), 2046--2070. — Dimensionally-split transport.
- Colella, P. and P. R. Woodward (1984). "The Piecewise Parabolic Method (PPM) for gas-dynamical simulations." *Journal of Computational Physics*, 54(1), 174--201. — PPM reconstruction and monotonicity limiter.
- Lin, S.-J. (2004). "A 'Vertically Lagrangian' Finite-Volume Dynamical Core for Global Models." *Monthly Weather Review*, 132(10), 2293--2307. — Vertically Lagrangian FV framework.
- Putman, W. M. and S.-J. Lin (2007). "Finite-volume transport on various cubed-sphere grids." *Journal of Computational Physics*, 227(1), 55--78. — FV transport on cubed-sphere grids.
- Ronchi, C., R. Iacono, and P. S. Paolucci (1996). "The 'Cubed Sphere': A New Method for the Solution of Partial Differential Equations in Spherical Geometry." *Journal of Computational Physics*, 124(1), 93--114. — Gnomonic cubed-sphere projection and metric tensors.
