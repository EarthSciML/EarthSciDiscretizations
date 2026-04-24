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
