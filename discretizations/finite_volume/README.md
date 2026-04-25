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
