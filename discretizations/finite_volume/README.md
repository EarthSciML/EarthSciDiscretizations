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
