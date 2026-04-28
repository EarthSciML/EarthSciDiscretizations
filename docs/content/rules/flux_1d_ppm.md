---
title: "flux_1d_ppm"
slug: "flux_1d_ppm"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(h^3) interior; O(h) at active limiter constraints"
applies_to: "flux(q, c, v) [form: ppm]"
rule_path: "discretizations/finite_volume/flux_1d_ppm.json"
description: "PPM-based 1D flux-form transport face flux: 4th-order CW84 edge interpolation + Colella-Woodward (1984) §4 monotonicity limiter + Courant-fraction flux integral + ifelse upwind selection on the per-face Courant."
tags: ["finite-volume", "ppm", "transport", "colella-woodward", "courant-integral"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/flux_1d_ppm-stencil.png"
       alt="6-point cartesian stencil at offsets {-3, -2, -1, 0, 1, 2} relative to the right cell of the interface F_{i+1/2}">
  <figcaption>6-point cartesian stencil at offsets <code>{-3, -2, -1, 0, 1, 2}</code>
  relative to the right cell of the interface <code>F_{i+1/2}</code>.
  Cells <code>{-2, -1, 0, 1}</code> drive the 4th-order CW84 edge interpolation
  <code>q_{i+1/2}</code>; the additional <code>{-3, +2}</code> cells are needed
  for the left-cell and right-cell parabolic reconstructions that feed the
  CW84 limiter and the Courant-fraction flux integral.</figcaption>
</figure>

## Continuous form

`flux_1d_ppm` discretizes the conservative 1D advective flux on a uniform
Cartesian grid,

$$F(x_{i+1/2},\,t)\;=\;v(x_{i+1/2},\,t)\,q(x_{i+1/2},\,t),$$

using the **piecewise-parabolic method** of Colella & Woodward (1984): each
cell carries a sub-grid parabolic profile reconstructed from cell-averaged
data, the §4 monotonicity limiter clips the parabola when an extremum sits
on the wrong side of the cell, and the face flux is the area-averaged
value of the upwind cell's profile over the swept volume `|c|·dx` of one
time step:

$$F_{i+1/2}\;=\;v_{i+1/2}\cdot
\begin{cases}
\frac{1}{|c|}\displaystyle\int_{1-|c|}^{1} a_{\text{left}}(\xi)\,d\xi & \text{if } c \ge 0,\\[1ex]
\frac{1}{|c|}\displaystyle\int_{0}^{|c|} a_{\text{right}}(\xi)\,d\xi & \text{if } c < 0.
\end{cases}$$

Upwind cell selection is encoded directly in the AST via
`ifelse(c >= 0, int_left, int_right)` — no caller-side branching on the
sign of the Courant number.

## Discrete operator

The closed-form composition is recorded in the JSON rule's
`edge_value_stencil`, `limiter`, `flux_integral`, and `flux_form` blocks.
At a high level:

1. **4th-order edge interpolation** (CW84 eq. 1.6) — same coefficients as
   [`ppm_reconstruction`]({{< ref "/rules/ppm_reconstruction" >}})'s
   `q_right_edge` stencil:
   $$q_{i+1/2}\;=\;\tfrac{7}{12}(q_{i-1}+q_i)\;-\;\tfrac{1}{12}(q_{i-2}+q_{i+1}).$$
   The same formula shifted by `±1` cell yields the left and right edges
   of each upwind cell's parabola (`ql_L_raw`, `qr_L_raw`, `ql_R_raw`,
   `qr_R_raw`).

2. **Colella & Woodward (1984) §4 monotonicity limiter** — applied
   independently to each upwind cell's reconstruction. The `is_extremum`
   case flattens to a constant when `q_i` is a local extremum within the
   stencil; the `overshoot_left` / `overshoot_right` cases pull the
   offending edge inward when the parabola's interior extremum sits past
   the opposite edge. Encoded as a closed-form `ifelse` AST identical to
   the [`vertical_remap`]({{< ref "/rules/vertical_remap" >}}) `limiter`
   block and bit-for-bit equivalent to
   [`_ppm_limit_cw84_sym`]({{< param repoURL >}}/blob/main/src/operators/reconstruction.jl).

3. **Courant-fraction flux integral** — closed-form parabolic
   antiderivative evaluated over the swept volume of the upwind cell.
   With `dq = qr - ql` and `q_6 = 6·(q_i − (ql + qr)/2)`,
   $$\text{int\_left}\;=\;qr - \tfrac{|c|}{2}\!\left[dq - q_6\!\left(1 - \tfrac{2}{3}|c|\right)\right],$$
   $$\text{int\_right}\;=\;ql + \tfrac{|c|}{2}\!\left[dq + q_6\!\left(1 - \tfrac{2}{3}|c|\right)\right].$$

4. **Upwind selection + face velocity** —
   $$F_{i+1/2}\;=\;v\cdot\text{ifelse}(c \ge 0,\;\text{int\_left},\;\text{int\_right}).$$

The full machine-readable AST lives in
[the JSON rule file]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_1d_ppm.json).

### Bindings

| Binding | Resolved by |
|---|---|
| `$q` | cell-averaged scalar field, ghost-extended by `Ng = 3` along the stencil axis |
| `$c` | per-face Courant number `c_{i+1/2} = u_{i+1/2}·dt/dx` (uniform Cartesian) or via the cubed-sphere panel-boundary-aware accessor (`_get_courant_xi/_get_courant_eta`) |
| `$v` | per-face velocity `v_{i+1/2}` |

The face-staggered `$c` and `$v` follow the same per-face binding contract
as [`lax_friedrichs_flux`]({{< ref "/rules/lax_friedrichs_flux" >}}) — the
analog of the slope-ratio `$r` in the limiter rules. The 6-point stencil
reaches three cells past the right cell of the interface; the imperative
reference in
[`src/operators/flux_1d.jl`]({{< param repoURL >}}/blob/main/src/operators/flux_1d.jl)
resolves this with a ghost-extended `q_ext`. ESS does not yet declare how
a rule consumes ghost cells — see the rule's `schema_gaps` block.

## Composition

`flux_1d_ppm` is the high-order PPM sibling of
[`lax_friedrichs_flux`]({{< ref "/rules/lax_friedrichs_flux" >}}): same
`op = flux` shape and per-face binding contract, but composes
[`ppm_reconstruction`]({{< ref "/rules/ppm_reconstruction" >}}) edge
values, the Colella-Woodward §4 monotonicity limiter (same closed-form
AST as the [`vertical_remap`]({{< ref "/rules/vertical_remap" >}})
`limiter` block), and a Courant-fraction flux integral, with `ifelse`
upwind selection in place of the LF `abs`-based first-order upwind. The
cell-centred tendency follows from the FV divergence
$T_i = -(F_{i+1/2}\cdot\Delta x_{i+1/2} - F_{i-1/2}\cdot\Delta x_{i-1/2})/A_i$
(or composes with
[`divergence_arakawa_c`]({{< ref "/rules/divergence_arakawa_c" >}}) on a
C-grid).

## Convergence

Numeric coverage today lives in the canonical Julia tests
[`test/test_transport_1d.jl`]({{< param repoURL >}}/blob/main/test/test_transport_1d.jl)
and
[`test/test_transport_2d.jl`]({{< param repoURL >}}/blob/main/test/test_transport_2d.jl)
(constant-field, linearity, conservation, and cubed-sphere advection
checks for `flux_1d_ppm!` and `flux_1d_ppm_arrayop`), and the
DCMIP-style cubed-sphere advection integration case in
[`test/integration_cases/cubed_sphere_advection.jl`]({{< param repoURL >}}/blob/main/test/integration_cases/cubed_sphere_advection.jl).

The hand-pinned single-face fixture at
[`fixtures/canonical/`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_1d_ppm/fixtures/canonical)
documents the closed-form face flux on a smooth sinusoidal profile
(positive- and negative-Courant cases; CW84 limiter inactive); the same
values are exercised inline in
[`test/test_flux_1d_ppm_rule.jl`]({{< param repoURL >}}/blob/main/test/test_flux_1d_ppm_rule.jl)
against `_ppm_limit_cw84_sym` + the closed-form Courant integral.

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The Layer-B walker harness needs (1) per-face Courant/velocity bindings
(<code>$c</code>, <code>$v</code>), (2) a ghost-extended input contract
(the 6-point stencil reaches three cells past the right cell of the
face), and (3) a Layer-B′ MMS-transport fixture kind — PPM flux + FV
divergence + RK time stepping — to verify convergence on a smooth
profile. Until those land, the convergence fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_1d_ppm/fixtures/convergence)
declares <code>applicable: false</code> and the rendered convergence
plot is suppressed.
</div>

## Reference

- Imperative reference: `flux_1d_ppm!` and `flux_1d_ppm_arrayop` in
  [`src/operators/flux_1d.jl`]({{< param repoURL >}}/blob/main/src/operators/flux_1d.jl);
  symbolic-tracing-safe limiter `_ppm_limit_cw84_sym` in
  [`src/operators/reconstruction.jl`]({{< param repoURL >}}/blob/main/src/operators/reconstruction.jl).
- Theory: Colella & Woodward (1984), JCP 54(1):174-201, eqs. (1.5)-(1.10)
  for the parabolic reconstruction and the (1.7)-(1.10) monotonicity
  limiter; the §4 Courant-fraction flux integral is the standard PPM
  upwind flux. Flux-form 1D advection composition follows Lin & Rood
  (1996) MWR.
- First-order sibling: [`lax_friedrichs_flux`]({{< ref "/rules/lax_friedrichs_flux" >}})
  (`F_{i+1/2} = max(c,0)·q_i + min(c,0)·q_{i+1}`).
- 2D cubed-sphere PPM transport sibling pending separate beads
  (catalog rows `fv3_lin_rood_advection` and
  `cam5_fv_ppm_reconstruction`).
