# Finite-Volume Rules

Rule files for finite-volume discretizations (PPM, MUSCL, WENO
reconstructions; flux forms; Riemann solvers, etc.).

The convention for file naming and organization within this directory may
evolve as content lands. Current expectation: one JSON file per named scheme,
e.g. `ppm_reconstruction.json`, `muscl_minmod.json`.

## Rules

- [`ppm_reconstruction.json`](ppm_reconstruction.json) — Colella & Woodward
  (1984) piecewise parabolic reconstruction, 1D, uniform Cartesian grid,
  unlimited (no monotonicity fix). Stencil i−2..i+2; 4th-order edge
  interpolation (eq. 1.6); parabola coefficients per eqs. (1.5),(1.7),(1.10).
  Layer-B MMS convergence at `tests/fixtures/ppm_reconstruction/` reaches
  ≥3.0 asymptotic order.
- [`weno5_advection.json`](weno5_advection.json) — Jiang & Shu (1996)
  classical WENO-5 flux reconstruction, 1D, uniform Cartesian grid,
  upwind-biased (caller picks `q^L` or `q^R` by local velocity sign).
  Stencil i−2..i+2 with three 3-point candidate sub-stencils (Shu 1998
  eq. 2.11), optimal linear weights d₀=1/10, d₁=6/10, d₂=3/10 (eq. 2.15),
  Jiang-Shu smoothness indicators β_k (eq. 2.17) and nonlinear ω_k (eqs.
  2.9–2.10). The ε-regularised ratio form of the nonlinear weights is not
  yet encodable in the ESS §7 stencil schema (follow-up ESS schema
  extension bead pending); the rule captures the symbolic form in a
  `nonlinear_weights` block. Layer-B fixtures at
  `tests/fixtures/weno5_advection/` exercise smooth MMS convergence
  (≥4.7 asymptotic order on phase-shifted sine) and a linear-advection
  square-wave shock-capturing sanity check (max overshoot/undershoot
  below 0.05 after one period).
- [`flux_limiter_minmod.json`](flux_limiter_minmod.json) — Roe (1986)
  minmod limiter, φ(r) = max(0, min(r, 1)). TVD, monotonicity-preserving,
  symmetric. The limiter **formula is encoded directly as an
  ExpressionNode AST** in the rule's `formula` field — there is no
  runtime-callable Julia name. Evaluation walks the AST. Applies to the
  slope-ratio scalar `$r`, so the rule is grid-family-agnostic in practice
  even though the POC scopes `grid_family` to `cartesian`. Layer-B
  fixture at `tests/fixtures/flux_limiter_minmod/phi_sweep.esm` pins the
  expected φ(r) at seven checkpoints and verifies Sweby-region bounds;
  `tvd_check.esm` confirms strict TVD on a smooth+square-wave IC over one
  advection period.
- [`flux_limiter_superbee.json`](flux_limiter_superbee.json) — Roe (1986)
  superbee limiter, φ(r) = max(0, min(2r, 1), min(r, 2)). Same shape and
  AST-first encoding as minmod; superbee sits on the upper edge of the
  Sweby (1984) second-order TVD region and is compressive near
  discontinuities. Layer-B fixtures parallel to minmod at
  `tests/fixtures/flux_limiter_superbee/`.
- [`transport_2d.json`](transport_2d.json) — first-order Rusanov / local
  Lax-Friedrichs flux-form 2D transport on the gnomonic cubed sphere
  with C-grid Courant bindings. 5-point in-panel stencil over offsets
  `(0, 0)`, `(±1, 0)`, `(0, ±1)`; coefficients depend on per-face
  Courant numbers (`c_xi_E`, `c_xi_W`, `c_eta_N`, `c_eta_S`), face
  lengths (`dx_E`, `dx_W`, `dy_N`, `dy_S`), and cell area (`A`) bound at
  the central cell by the cubed_sphere accessor — the same
  accessor-resolves-staggered-bindings convention as
  [`covariant_laplacian_cubed_sphere`](../finite_difference/covariant_laplacian_cubed_sphere.json).
  Mirrors the imperative `transport_2d` ArrayOp in
  `src/operators/transport_2d.jl`. Layer-A canonical fixture under
  `transport_2d/fixtures/canonical/` documents the bit-equivalence
  contract against `transport_2d(q, courant_xi, courant_eta,
  CubedSphereGrid(24))` to within `1e-12`. Layer-B convergence declares
  `applicable: false` — ESS's `verify_mms_convergence` cannot drive a 2D
  cubed-sphere staggered C-grid sweep with face-Courant bindings today;
  numeric coverage continues to live at `test/test_transport_2d.jl` and
  `test/integration_cases/cubed_sphere_advection.jl`.
- [`divergence_arakawa_c.json`](divergence_arakawa_c.json) — Arakawa &
  Lamb (1977) C-grid divergence: ∂uₓ/∂x + ∂u_y/∂y with uₓ at face-x and
  u_y at face-y, output at cell center. First rule to use the
  `arakawa` selector kind (`stagger` ∈ {`cell_center`, `face_x`,
  `face_y`, `vertex`}) and the file-local `enums` block (ESS §9.3) to
  carry stagger values portably across bindings. Two-point centered
  per axis, O(h²). The schema decisions for the new `kind` and the
  stagger-enum convention are pinned in
  [`../SELECTOR_KINDS.md`](../SELECTOR_KINDS.md) (rows 6–9).
  Layer-A fixture under `divergence_arakawa_c/fixtures/canonical/`
  pins a 2×2 cell example (u(x,y) = x at face_x, v(x,y) = y² at
  face_y → div = [[1.5, 2.5], [1.5, 2.5]]). Layer-B convergence
  declares `applicable: false` — ESS's `verify_mms_convergence` is
  wired only for 1D periodic Cartesian stencils today, so the 2D
  staggered MMS sweep awaits an ESS harness extension; numeric
  coverage continues to live alongside the runtime tests.
- [`vertical_remap.json`](vertical_remap.json) — Lin (2004) MWR
  PPM-based conservative vertical remap with Colella-Woodward (1984)
  monotonicity limiting. Each old layer's parabolic profile is
  reconstructed (interior 4-point CW84 stencil; 1st/2nd-order near
  the column top and bottom), CW84-limited (encoded as a closed-form
  ifelse-based AST under the rule's `limiter` block), and the
  closed-form parabolic antiderivative is summed across the
  cumulative-pressure intersection of the new layers with the old
  layers. The rule is *conservative* (column-integrated `q · dp`
  is preserved exactly in real arithmetic) and *monotonicity-
  preserving* under the CW84 limiter.

  **Status: documentation-only — deferred to a future phase-hooks
  RFC (dsc-otd).** Conservative remap is structurally a phase-hook
  operation (Lagrangian → Eulerian re-gridding between timesteps)
  with a variable-arity, data-dependent neighbor list and time-
  varying per-column metrics, none of which fit the §7 stencil
  rule schema. ESS will introduce a `phase_hooks` /
  `integration_callbacks` schema block in a future RFC; until then
  the rule file is retained as a reference artifact (math,
  references, AST fragments) and is *not* executed by the rule
  engine. There is **no imperative implementation in any binding
  (Julia/Python/Rust/TypeScript) and none is to be re-added** —
  see the rule's `schema_gaps` block for the disposition. The
  walker discovers the rule (Layer-A canonical-form discovery
  exercises the file) and reports Layer-B as
  "fixture-declared not applicable" with the phase-hook deferral
  reason (see [`Lagrangian vertical scheme`](#lagrangian-vertical-scheme--vertical-remapping-not-yet-executable)
  below).
- [`lax_friedrichs_flux.json`](lax_friedrichs_flux.json) — Lax &
  Friedrichs (1954) numerical flux for linear advection. Two-point
  cartesian flux stencil at the face: F_{i+1/2} = max(c,0)·q_i +
  min(c,0)·q_{i+1}, with the face-staggered Courant `$c` carried as a
  per-face binding (analogous to `$r` in the limiter rules) and upwind
  selection encoded directly in the AST via the `abs` op — no caller-side
  branching on sign(c). Reduces to first-order upwinding; dissipative by
  construction and retained as a debugging / oracle scheme. The cubed-
  sphere wrapper in `src/operators/flux_1d.jl` matches the in-panel core
  encoded here verbatim; cross-panel ghost extension and panel-boundary
  distance handling await schema follow-ups (boundary_policy +
  time-varying / face-stagger bindings, tracked off dsc-35x). Layer-A
  canonical fixture is intentionally omitted because ESS's `discretize`
  does not yet support op="flux" with per-face bindings — the
  hand-pinned 5-interior-face example (q = powers of 2; mixed-sign
  Courant including c=0) lives inline in
  `test/test_lax_friedrichs_flux_rule.jl` instead. Layer-B convergence
  declares `applicable: false` for the same face-stagger /
  per-face-binding gap as the limiter rules — a follow-up Layer-B′
  MMS-transport fixture (LF + divergence + forward Euler) is the right
  shape.

## Composing a limiter with a reconstruction

The limiter rules accept a slope-ratio scalar `$r` and return φ(`$r`).
They do not carry stencils — the slope ratio is the caller's
responsibility.

Worked example (1D uniform Cartesian, u > 0):

1. **Compute the slope ratio at each cell interface**:

   `r_i = (q_i - q_{i-1}) / (q_{i+1} - q_i + ε)` (guard against zero
   denominators with a small ε when the solution is locally flat).

2. **Evaluate the limiter rule's `formula` AST** at `$r = r_i` using the
   ESD evaluator (`EarthSciDiscretizations.eval_coeff`) or an equivalent
   tree-walk evaluator in your binding. The result is φ(r_i).

3. **Scale the high-order slope correction by φ(r_i)** in the MUSCL-style
   flux, e.g.

   `F_{i+1/2} = u * (q_i + 0.5 * φ(r_i) * (q_{i+1} - q_i))`.

   Under forward Euler with CFL ≤ 1/(1 + 0.5·φ_max) (φ_max = 1 for
   minmod, 2 for superbee), this scheme is strictly TVD. The Layer-B′
   fixtures in `tests/fixtures/flux_limiter_{minmod,superbee}/tvd_check.esm`
   drive exactly this composition with a 64-cell grid and CFL=0.4 to
   verify total variation is non-increasing over one advection period.

Combining a limiter with PPM (`ppm_reconstruction.json`) or WENO-5
(`weno5_advection.json`) follows the same pattern: compute `r` at each
interface from the relevant reconstruction values and scale the
high-order contribution by φ(`r`). The limiter rule is entirely
reconstruction-agnostic.

## Lagrangian vertical scheme — vertical remapping not yet executable

Atmospheric finite-volume dynamical cores using a Lagrangian vertical
coordinate (FV3-style) require periodic conservative remapping of fields
from the drifted Lagrangian layers back to a target Eulerian (hybrid
sigma-pressure) vertical grid. Conservative remap is structurally a
phase-hook operation — it runs between timesteps, has a data-dependent
variable-arity neighbor list (cumulative-pressure overlap of `dp_old`
with `dp_new`), and does not fit the §7 stencil rule schema. ESD's
current rule engine and walker operate without this step.

The math is documented declaratively in
[`vertical_remap.json`](vertical_remap.json) as a reference artifact
(the same way `weno5_advection` was carried before AST-walker dispatch
landed in ESS). The rule file is **not** executed today and there is
**no imperative implementation in any binding**. Models needing a
Lagrangian vertical configuration must either avoid the Lagrangian
path or use the Eulerian-vertical configurations supported today
(e.g. [`centered_2nd_uniform_vertical`](../finite_difference/centered_2nd_uniform_vertical.json)).

Tracking: a future ESS RFC will introduce a `phase_hooks` /
`integration_callbacks` schema block for declarative phase-hook
contracts (vertical remap, coupler exchange, restart IO). Until then
vertical remap is intentionally documentation-only. See
[`vertical_remap.json`](vertical_remap.json)'s `schema_gaps` block
and ESD/dsc-otd for the disposition.
