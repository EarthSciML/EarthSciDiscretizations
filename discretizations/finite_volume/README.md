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
  classical WENO-5 reconstruction, 1D, uniform Cartesian grid,
  upwind-biased. **Canonical nonlinear-scheme exemplar (dsc-b78):**
  `applies_to` matches the §4.2 op `div(*($U, $q))` and the lowering is
  a single closed `arrayop` whose body is the FV divergence
  `(F_E − F_W)/dx` with `F = U_face · ifelse(U_face > 0, q^L, q^R)`.
  The Jiang-Shu smoothness indicators β_k (eq. 2.17), ε-regularised
  nonlinear weights α_k = d_k/(ε+β_k)² (eqs. 2.9–2.10), and the convex
  combination ω_k·p_k are all expressed as `+`, `-`, `*`, `/`, `^` over
  `index` selectors — no scheme-specific `stencil` / `selector` /
  `offset` blobs and no off-spec match keys. The right-biased
  reconstruction `q^R` is the same closed expression evaluated on the
  mirrored 5-cell stencil (Shu 1998 §2.2 eq. 2.16); ESS evaluates the
  AST in any binding by walking the existing arrayop / broadcast
  machinery. Layer-B fixtures at `tests/fixtures/weno5_advection/`
  exercise smooth MMS convergence (≥4.7 asymptotic order on
  phase-shifted sine) and a linear-advection square-wave
  shock-capturing sanity check (max overshoot/undershoot below 0.05
  after one period).
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
- [`ppm_edge_cubed_sphere.json`](ppm_edge_cubed_sphere.json) — Harris et
  al. (2021), GFDL FV3 Technical Memorandum, Eq. 6.5–6.6: two-sided PPM
  extrapolation at cubed-sphere panel-boundary interfaces. The standard
  4th-order PPM edge formula breaks down across the discontinuous panel
  coordinate system; this rule averages two one-sided extrapolations
  using only in-panel offsets and lets the cubed_sphere grid accessor
  (`src/grids/panel_connectivity.jl`) resolve cross-panel offsets to
  neighbor-panel ghost cells. Single-axis stencil at offsets
  `{-1, 0, 1, 2}` along `$x ∈ {xi, eta}`; coefficients
  `[-1/4, 3/4, 3/4, -1/4]` after specializing the FV3 non-uniform
  formula to the gnomonic equidistant grid's uniform isotropic spacing
  `h = π/(2·Nc)`. The FV3 eq. 6.6 monotonicity clamp is documented in
  the rule's `monotonicity_constraint` block as a non-linear
  post-processing step (post-stencil clamps are a follow-up ESS schema
  extension; see flux_limiter_minmod's `formula` AST for the related
  pattern). Layer-A canonical fixture pins bit-equivalence with the
  pre-port `ppm_edge_value_twosided` reference on `c24`; Layer-B
  convergence declares `applicable: false` until ESS gains
  cubed_sphere selectors with panel-connectivity dispatch (shared
  follow-up with `covariant_laplacian_cubed_sphere`).
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
  preserving* under the CW84 limiter. The conservative-remap rule
  type and the time-varying per-cell vertical metric binding are two
  ESS schema extensions that this rule needs before it can be
  materialized straight from JSON; both are documented in the rule's
  `schema_gaps` block (see also `docs/rule-catalog.md` row
  `vertical_remap_lin_2004` and `discretizations/SELECTOR_KINDS.md`
  decision #3). Walker fixtures under
  `vertical_remap/fixtures/{convergence,conservation,monotonicity}/`
  declare `applicable: false` and link to the runtime coverage at
  `test/test_vertical_remap.jl` until the schema gaps land. The
  Julia runtime path at `src/operators/vertical_remap.jl` remains
  the live implementation.
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
- [`lax_friedrichs_flux_cubed_sphere_xi.json`](lax_friedrichs_flux_cubed_sphere_xi.json)
  / [`lax_friedrichs_flux_cubed_sphere_eta.json`](lax_friedrichs_flux_cubed_sphere_eta.json)
  — sibling-pair cubed-sphere face-stagger wrappers for the LF flux,
  one rule per axis (ξ → emits at u_edge, η → emits at v_edge).
  Two-point in-panel stencil at offsets `{-1, 0}` along the chosen
  axis at `cell_center` stagger; the face-staggered Courant `$c` is
  consumed via a `reads` block at u_edge / v_edge offset 0, matching
  the pattern from `fv3_sinsg_flux_xi.json`. The closed-form
  upwind-blended F = (c+|c|)/2·q_west + (c-|c|)/2·q_east mirrors the
  Cartesian core (`lax_friedrichs_flux.json`) coefficient-for-
  coefficient — only selector kind / stagger differs. Cross-panel
  ghost extension and panel-boundary distance handling for the
  Courant precomputation live in the cubed_sphere grid accessor
  (`src/grids/panel_connectivity.jl` + `_get_courant_xi/eta` in
  `src/operators/flux_1d.jl`); selectors carry no `panel` field per
  SELECTOR_KINDS.md decision #13. Layer-A canonical fixture pins
  per-face flux equivalence on c4 against the imperative `flux_1d`
  reference within `1e-12`; Layer-B convergence is `applicable: false`
  pending the cubed_sphere walker dispatch + face-stagger MMS-transport
  harness; Layer-C integration carries a Williamson-1 cosine-bell stub
  behind ESD_RUN_INTEGRATION=1. The η rule's Layer-A fixture is
  skip-only and points at the ξ fixture (the algebra is symmetric
  under coordinate reflection). Tracked at dsc-0fd; the `boundary_policy`
  schema feature that would let these rules reference
  `grid.dist_xi_bnd` / `grid.dist_eta_bnd` directly is a follow-up.

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
