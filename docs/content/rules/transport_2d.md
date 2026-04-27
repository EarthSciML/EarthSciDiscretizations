---
title: "transport_2d"
slug: "transport_2d"
families: "finite_volume"
grid_families: "cubed_sphere"
rule_kinds: "scheme"
accuracy: "O(h)"
applies_to: "advect(q, U)"
rule_path: "discretizations/finite_volume/transport_2d.json"
description: "First-order Rusanov / local Lax-Friedrichs flux-form 2D transport on the gnomonic cubed sphere with C-grid Courant bindings."
tags: ["finite-volume", "cubed-sphere", "transport", "lax-friedrichs", "rusanov"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/transport_2d-stencil.png"
       alt="5-point in-panel stencil with face-Courant-weighted coefficients">
  <figcaption>5-point in-panel stencil over the <code>(ξ, η)</code> coordinates.
  The center carries the inflow contribution from every face; each cardinal
  neighbor carries the inflow / outflow contribution from a single face,
  weighted by the local Courant number and physical face length.</figcaption>
</figure>

## Continuous form

`transport_2d` discretizes the conservative 2D advection operator on a
panel of the gnomonic cubed sphere,

$$\frac{\partial q}{\partial t} \;=\; -\frac{1}{A}\bigl[\partial_\xi(u\,q\,\Delta x) + \partial_\eta(v\,q\,\Delta y)\bigr],$$

using a first-order **local Lax-Friedrichs (Rusanov)** flux at each face
where the dissipation coefficient is set to the local face Courant
magnitude:

$$F_{\text{face}} \;=\; \frac{c}{2}\,(q_L + q_R) \;-\; \frac{|c|}{2}\,(q_R - q_L).$$

With `α = |c|` the LF flux reduces to first-order donor-cell upwind, so
this rule is equivalent to first-order upwind on every interface but is
named `transport_2d` to match the imperative `transport_2d` ArrayOp
retained in [`src/operators/transport_2d.jl`]({{< param repoURL >}}/blob/main/src/operators/transport_2d.jl)
for symbolic / debugging use alongside the higher-order Lin-Rood and
unsplit-PPM transport paths.

## Discrete operator

Expanding the symmetric and dissipative parts of the LF flux gives a
linear combination of the central cell and its four cardinal neighbors:

$$\frac{dq_{i,j}}{dt} \;=\; \sum_{(\Delta\xi,\Delta\eta)\in\{(0,0),\,(\pm1,0),\,(0,\pm1)\}} c_{\Delta\xi,\Delta\eta}\;q_{i+\Delta\xi,\,j+\Delta\eta},$$

with face-Courant-dependent coefficients (cell area `A`, face lengths
`dx_E`, `dx_W`, `dy_N`, `dy_S`):

| offset `(Δξ, Δη)` | coefficient `c_{Δξ,Δη}` |
|---|---|
| `( 0,  0)` | `−[((c_E+|c_E|)/2)·dx_E − ((c_W−|c_W|)/2)·dx_W + ((c_N+|c_N|)/2)·dy_N − ((c_S−|c_S|)/2)·dy_S] / A` |
| `(+1,  0)` | `−((c_E−|c_E|)/2)·dx_E / A` (active when `c_E < 0`) |
| `(−1,  0)` | `+((c_W+|c_W|)/2)·dx_W / A` (active when `c_W > 0`) |
| `( 0, +1)` | `−((c_N−|c_N|)/2)·dy_N / A` (active when `c_N < 0`) |
| `( 0, −1)` | `+((c_S+|c_S|)/2)·dy_S / A` (active when `c_S > 0`) |

The full machine-readable coefficient AST lives in
[the JSON rule file]({{< param repoURL >}}/blob/main/discretizations/finite_volume/transport_2d.json).

### Bindings

| Binding | Resolved by |
|---|---|
| `c_xi_E`, `c_xi_W` | cubed_sphere C-grid accessor at the east / west faces of the central cell |
| `c_eta_N`, `c_eta_S` | cubed_sphere C-grid accessor at the north / south faces of the central cell |
| `dx_E`, `dx_W`, `dy_N`, `dy_S` | physical face lengths from `grid.dx`, `grid.dy` at the four faces of the central cell |
| `A` | cell area from `grid.area` at the central cell |

Cross-panel ghost cells and basis rotation at panel boundaries live in
the cubed_sphere accessor (`src/grids/panel_connectivity.jl`); the rule
selectors do **not** carry a `panel` field — `(Δξ, Δη)` offsets resolve
across panel seams via the accessor's connectivity table. This is the
same convention as
[`covariant_laplacian_cubed_sphere`]({{< ref "/rules/covariant_laplacian_cubed_sphere" >}}).

## Convergence

Numeric coverage today lives in the canonical Julia tests
[`test/test_transport_2d.jl`]({{< param repoURL >}}/blob/main/test/test_transport_2d.jl)
(constant-field, linearity, and steady-state checks for the
`transport_2d` ArrayOp) and the cubed-sphere advection integration case
in
[`test/integration_cases/cubed_sphere_advection.jl`]({{< param repoURL >}}/blob/main/test/integration_cases/cubed_sphere_advection.jl).
The bit-equivalence canonical fixture at
[`fixtures/canonical/`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/transport_2d/fixtures/canonical)
documents the equivalence contract against `transport_2d(q,
courant_xi, courant_eta, CubedSphereGrid(24))` to within `1e-12`
relative tolerance.

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The Layer-B walker harness needs the in-flight 2D dispatch + per-cell
C-grid face-staggered callables (Courant numbers + face lengths) to
evaluate this rule on a manufactured solution defined in
<code>(ξ, η)</code>. Until that extension lands, the convergence
fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/transport_2d/fixtures/convergence)
declares <code>applicable: false</code> and the rendered convergence
plot is suppressed.
</div>

## Reference

- Imperative reference: `transport_2d` in
  [`src/operators/transport_2d.jl`]({{< param repoURL >}}/blob/main/src/operators/transport_2d.jl).
- Higher-order siblings on the same grid: Lin-Rood (`transport_2d_linrood!`)
  and unsplit-PPM (`transport_2d_ppm_arrayop`) — both pending separate
  beads (catalog rows `fv3_lin_rood_advection` and
  `cam5_fv_ppm_reconstruction`).
