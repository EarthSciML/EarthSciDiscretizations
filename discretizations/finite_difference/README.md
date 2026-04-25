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
  (uniform sigma; non-uniform vertical metrics are deferred — see
  `SELECTOR_KINDS.md` decision #3).
- `centered_2nd_uniform_latlon.json` — latlon, axes `lon` / `lat` (literal
  per `SELECTOR_KINDS.md` decision #6), angular spacings `dlon` / `dlat`,
  sphere radius `R`, latitude metric `cos_lat`. The lon-axis coefficient
  carries `1/(R cos_lat dlon)`; the lat-axis coefficient carries `1/(R dlat)`.
  Layer B marks `applicable: false` pending ESS MMS-harness support for 2D
  structured stencils + per-cell metric bindings (decision #8 + follow-up).
