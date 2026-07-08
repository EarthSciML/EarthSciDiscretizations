# Finite-Difference Rules

Rule files for finite-difference discretizations on structured grids
(centered, upwind, compact stencils, etc.).

The convention for file naming and organization within this directory may
evolve as content lands. Current expectation: one JSON file per named scheme,
e.g. `central_4th_order.json`, `upwind_3rd_order.json`.

## Naming across grid families

When the same scheme exists for multiple grid families (cartesian, vertical,
lat-lon, …), a per-family suffix is appended to the rule name and a sibling
file is authored. The selector `kind` field inside the rule names the grid
family explicitly. See `../SELECTOR_KINDS.md` for the running design index.

Current variants:

- `centered_2nd_uniform.json` — cartesian, axis `$x`, spacing `dx`/`h`.
- `centered_2nd_uniform_vertical.json` — vertical, axis `$k`, spacing `h`
  (uniform sigma).
- `centered_2nd_nonuniform_vertical.json` — vertical, axis `$k`, per-cell
  layer thickness `dz[k]` via `{op: "index", args: ["dz", "$k"]}` — see
  `../SELECTOR_KINDS.md` decision #3.
- `centered_2nd_uniform_latlon.json` — latlon, axes `lon` / `lat` (literal
  per `SELECTOR_KINDS.md` decision #10), angular spacings `dlon` / `dlat`,
  sphere radius `R`, latitude metric `cos_lat`. The lon-axis coefficients
  carry `±1/(2 R cos_lat dlon)`; the lat-axis coefficients carry
  `±1/(2 R dlat)`. The Layer-B convergence fixture is active (Y_{2,0}
  spherical-harmonic MMS sweep, n = 16→128) and driven by the walker's
  `2d_latlon_sphere` runner (`test/walk_esd_tests.jl`).
- `covariant_laplacian_cubed_sphere.json` — cubed_sphere, axes `xi`/`eta`,
  spacing `h`. Encoded as a `kind: "cross_metric"` composite: 8 `terms`
  pair helper axis-stencils (`d2_dxi2_cubed_sphere`,
  `d2_deta2_cubed_sphere`, `d2_dxieta_cubed_sphere`,
  `d1_dxi_over_J_cubed_sphere`, `d1_deta_over_J_cubed_sphere` — defined
  in the same file) with per-cell metric components (`ginv_xi_xi`,
  `ginv_eta_eta`, `ginv_xi_eta`, `dJgxx_dxi`, `dJgyy_deta`, `dJgxe_dxi`,
  `dJgxe_deta`). The helper sub-schemes carry the per-entry `selectors`
  arrays (one per axis) and the `h` / `J` metric bindings. Cross-panel
  ghost cells are resolved by the cubed_sphere grid accessor — selectors
  do **not** carry a `panel` field. See `SELECTOR_KINDS.md` decision #13.
