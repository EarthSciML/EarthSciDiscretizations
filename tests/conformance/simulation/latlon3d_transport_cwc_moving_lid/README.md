# The moving-lid consistency-with-continuity gate

Free-stream preservation for the **terrain-following** vertical transport stack while
the level thicknesses **breathe in time** — the production invariant for a
hybrid sigma-pressure (GEOS-FP-style) vertical coordinate. It is the moving-grid
counterpart of `latlon3d_transport_cwc_3d_ppm`, which holds the vertical geometry
fixed.

## The finding behind it

The fixed-grid vertical flux divergence is `D = (M_{k+1/2} − M_{k−1/2}) / dz[k]`. That
`/dz` makes the prognostic a pseudo-density and is conservative **only when `dz` is
constant in time**: on a breathing grid it omits the grid-motion (thickness-tendency)
term, so a uniform tracer drifts. The continuous layer budget is mass-form,
`∂(dp)/∂t = −∇·(v dp) − (M_{k+1/2} − M_{k−1/2})`, with **no** division.

So a terrain-following scheme must **separate two roles** the fixed-grid rule conflates
into one `dz`:

- the **divergence measure** → mass form, the raw flux difference (no `/dz`), so the
  cell air mass is integrated as an extensive quantity and conserves while the layers
  move; and
- the **reconstruction width** → the physical, breathing pressure thickness
  `dp[i,j,k,t]` (rank-3, terrain-following), the Lin–Rood mass-coordinate reconstruction.

This case exercises the two new rules that do exactly that —
`face_flux_divergence_lev_massform_noflux_bc` (air mass) and
`ppm_flux_D_lev_mono_hybrid_noflux_bc` (tracer) — with the horizontal (lon) rules
unchanged, since a lat–lon mesh only moves vertically.

## What it asserts

The model carries the **air mass** `m` (the breathing pressure thickness) and the
**tracer mass** `mq`, derives `q = mq/m`, and reconstructs the vertical PPM on the
breathing widths `dp = m`. The vertical air-mass flux `Mz` is **diagnosed from
continuity** in closed form — `Mz[i,j,kn] = −horiz[i,j,1]·φ(kn)` with
`φ(kn) = (kn−1) − (kn−1)²/NLEV` vanishing at both walls — so the mass-form continuity
reproduces the hybrid tendency `dB[k]·∂Ps/∂t` exactly and `m` stays equal to
`dA[k] + dB[k]·Ps(i,j,t)` as the surface pressure redistributes zonally.

At `q == 1` every CW84 slope is exactly `0.0`, so the parabola is the constant 1
**regardless of the breathing widths**, the flux collapses to `F = Mz` bitwise, and the
terrain-following PPM collapses onto the mass-form continuity. So `dmq/dt` equals
`dm/dt` to the last bit and **`dev = mq − m` stays exactly 0 while the grid moves**.

**Measured** (julia, `Tsit5`, reltol 1e-10 / abstol 1e-12): `dev` is exactly `0.0` at
t=0 and t=0.05. The layers genuinely breathe — `m` moves by `max|Δ| ≈ 7.2e-3` (~1.5%)
over the span and stays strictly positive (min ≈ 0.449) — and the reconstruction is
confirmed terrain-following: the lowered right-hand side carries 25,441 rank-3 `dp`
width gathers and the diagnosed `Mz`.

## Two things to know

**1. The reference is the constant zero, so `0.0` proves nothing by itself.**
Non-vacuity is established by `problems/_falsifier_moving_lid.esm` — the same problem
with `ic(mq)` multiplied by 1.000001, so `q = 1+1e-6`. With `q ≠ 1` the parabola is no
longer flat, the tracer stops tracking the breathing air mass, and `dev` reaches
**7.24e-9** at t=0.05 — six orders above the 1e-12 tolerance — and the case correctly
**fails**. The drift is exactly `(breathing 7.24e-3) × (perturbation 1e-6)`, the
expected `(q−1)·(air-mass tendency)` rate. The falsifier is regenerated from this
problem plus the single `ic(mq)` node, so the two cannot drift apart; re-run it if you
change this problem.

**2. Only the t>0 assertion discriminates.** `dev` is a state integrated from
`dev(0) = 0`, so its t=0 value is `0.0` by construction even under perturbation.

## Cost and mesh

12×7×9 (NLEV=9 ≥ 7 for the monotone lev rule; NLON=12 for the zonal wrap; no lat
transport, so the mono-PPM box count stays small and the case builds in well under a
minute). The wind is a steady, cos(lat)-tapered zonal air-mass flux that drives a zonal
surface-pressure redistribution; the air mass breathes over the span while staying
positive.
