---
title: "weno5_advection"
slug: "weno5_advection"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx⁵)"
applies_to: "div(U·q), dim=x"
rule_path: "discretizations/finite_volume/weno5_advection.json"
description: "Jiang–Shu (1996) 5th-order WENO reconstruction expressed as a closed §4.2 arrayop+ifelse lowering — canonical nonlinear-scheme exemplar."
tags: ["finite-volume", "weno", "jiang-shu", "advection", "high-order"]
---

## Lowering

The rule is the canonical exemplar of an ESD **nonlinear** scheme that
lowers a §4.2 PDE operator to a **closed `arrayop` expression** in the
§4.2 op vocabulary — the sibling of
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}) for
linear schemes. There are no scheme-specific coefficient blobs and no
off-spec match keys: `applies_to` matches against `div` (a §4.2 op)
applied to the advective flux `*($U, $q)`, and the lowering is a single
`arrayop` whose body combines `index`, `+`, `-`, `*`, `/`, `^`, and
`ifelse`. ESS bindings (any binding) can evaluate the lowering by
walking the AST through the existing arrayop / broadcast evaluator —
no new kernels.

```text
applies_to:  div(*($U, $q), dim=x)
replacement: arrayop(
               output_idx = [$x],
               expr       = (F_E − F_W) / dx,
               args       = [$U, $q]
             )
```

The face flux `F = U_face · q^upwind` is built from a cell-to-face
average of `U` and a Jiang–Shu (1996) WENO5 reconstruction of `q`.
Upwinding is encoded directly in the AST as

$$F_{i+1/2} \;=\; U_{i+1/2}\;\cdot\; \mathrm{ifelse}\!\left(U_{i+1/2} > 0,\; q^{L}_{i+1/2},\; q^{R}_{i+1/2}\right),$$

with the cell-to-face velocity average $U_{i+1/2} = \tfrac12(U_i + U_{i+1})$
and a symmetric expression at the west face. The right-biased
reconstruction $q^{R}$ is the same closed WENO5 expression evaluated on
the mirrored 5-cell stencil (Shu 1998 §2.2 eq. 2.16) — the algebra is
identical, only the input order differs.

## Reconstruction (Jiang–Shu 1996)

For face $x_{i+1/2}$ with $U_{i+1/2} > 0$, the left-biased reconstruction
is the convex combination

$$q^{L}_{i+1/2} \;=\; \omega_0 p_0 + \omega_1 p_1 + \omega_2 p_2,$$

of three candidate sub-stencil polynomials evaluated at the face
(Shu 1998 eq. 2.11):

| sub-stencil | offsets | polynomial $p_k(q)$ |
|---|---|---|
| `p₀` | (−2, −1, 0)   | $\tfrac{1}{3}q_{i-2} - \tfrac{7}{6}q_{i-1} + \tfrac{11}{6}q_i$ |
| `p₁` | (−1, 0, +1)   | $-\tfrac{1}{6}q_{i-1} + \tfrac{5}{6}q_i + \tfrac{1}{3}q_{i+1}$ |
| `p₂` | (0, +1, +2)   | $\tfrac{1}{3}q_i + \tfrac{5}{6}q_{i+1} - \tfrac{1}{6}q_{i+2}$ |

The nonlinear weights $\omega_k$ are normalized smoothness-weighted
versions of the optimal linear weights
$(d_0, d_1, d_2) = (1/10, 6/10, 3/10)$:

$$\beta_0 = \tfrac{13}{12}(q_{i-2}-2q_{i-1}+q_i)^2 + \tfrac14(q_{i-2}-4q_{i-1}+3q_i)^2,$$
$$\beta_1 = \tfrac{13}{12}(q_{i-1}-2q_i+q_{i+1})^2 + \tfrac14(q_{i-1}-q_{i+1})^2,$$
$$\beta_2 = \tfrac{13}{12}(q_i-2q_{i+1}+q_{i+2})^2 + \tfrac14(3q_i-4q_{i+1}+q_{i+2})^2,$$
$$\alpha_k = \frac{d_k}{(\varepsilon + \beta_k)^2}, \qquad \omega_k = \frac{\alpha_k}{\sum_j \alpha_j}.$$

Every term above appears in the rule's `replacement` AST as `+`, `-`,
`*`, `/`, and `^` over `index` selectors — there is no
`smoothness_indicators` blob, no `nonlinear_weights` blob, and no
ratio-form ω that has to live outside the AST.

## Boundary conditions

Boundary handling is read from the domain's `boundary_conditions` block
(esm-spec §11.5) at lowering time and applied as **downstream rewrite
rules over concrete indices** — there is no `bc:*` op embedded in the
lowered AST. The pattern is:

| Domain BC | Index transformation applied to `$q[$x ± k]`, `$U[$x ± k]` |
|---|---|
| `periodic` | `mod($x ± k + N, N)` (see [`periodic_bc`]({{< ref "/rules/periodic_bc" >}})) |
| `dirichlet` / `constant` | Boundary cells read the prescribed value via `index` into a fill row |
| `zero_gradient` / `neumann` | Mirror the in-range neighbor (`min`/`max` clamp on the index) |

Each transformation is a separate rule that fires at the boundary
cells. The `weno5_advection` rule itself stays BC-agnostic — its
replacement is the interior closed form. The lowering pipeline (rule
application + BC rewrites) takes the `(grid_family, BC list)` pair from
the domain and emits index expressions that respect the declared BCs.

## Stencil

<figure class="figure">
  <img src="/plots/rules/weno5_advection-stencil.png"
       alt="WENO5 three sub-stencils, left-biased branch">
  <figcaption>Three candidate 3-cell sub-stencils for the left-biased
  reconstruction of the edge value <code>q<sub>i+1/2</sub><sup>L</sup></code>.
  Linear weights <code>(d₀, d₁, d₂) = (1/10, 6/10, 3/10)</code> recover formal
  5th-order accuracy in smooth regions; nonlinear weights ω<sub>k</sub> shift
  away from sub-stencils that straddle a discontinuity.</figcaption>
</figure>

For a single divergence at cell `i`, the closed AST reads from cells
`i−3, …, i+3` (seven cells) for `q` — five per face — and from cells
`i−1, i, i+1` for the cell-to-face velocity average. The right-biased
branch (used when `U_face < 0`) is the mirror of the left-biased branch
under index reflection.

## Convergence

Numeric coverage lives in the canonical Julia test fixture at
[<code>tests/fixtures/weno5_advection/</code>]({{< param repoURL >}}/blob/main/tests/fixtures/weno5_advection),
exercised by
[<code>test/test_weno5_advection_rule.jl</code>]({{< param repoURL >}}/blob/main/test/test_weno5_advection_rule.jl).
The MMS fixture uses `f(x) = sin(2πx + 1)` on `[0, 1]` with periodic boundary
conditions; the phase shift keeps the critical points of `f` away from every
dyadic cell face at `n ∈ {32, 64, 128, 256}`, sidestepping the well-known
WENO5-JS accuracy dip from `ω_k → d_k` recovery stalling at
`f'(x_{i+1/2}) = 0` (Henrick, Aslam & Powers, JCP 2005).

| `n`  | `dx`         | L∞ error     | observed order |
| ---: | :----------- | :----------- | :------------- |
|  32  | 0.03125000   | 3.4137e-05   | —              |
|  64  | 0.01562500   | 1.0656e-06   | 5.002          |
| 128  | 0.00781250   | 3.3254e-08   | 5.002          |
| 256  | 0.00390625   | 1.0377e-09   | 5.002          |

Theoretical asymptotic order: **5.0** (Jiang & Shu 1996, smooth regions).
Acceptance threshold: **min(observed order) ≥ 4.7** — leaves headroom for the
small accuracy hit from the `ε = 1e-6` regularisation of the nonlinear
weights. A companion shock-capturing fixture advects a unit square wave one
full period at CFL 0.4 with SSP-RK3; max overshoot/undershoot is ~3.7e-4,
well under the 0.05 tolerance.

The walker-side fixture under
[<code>discretizations/finite_volume/weno5_advection/fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/weno5_advection/fixtures/convergence)
exercises Layer-B convergence through ESS's AST-walker dispatch over the
closed `arrayop` + `ifelse` replacement — no scheme-specific kernels.
