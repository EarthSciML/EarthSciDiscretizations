# The regional lateral-inflow consistency-with-continuity gate

Free-stream preservation for the **regional (limited-area)** monotone-PPM horizontal
transport stack with **open lateral boundaries** carrying a prescribed Dirichlet inflow
concentration — the production invariant for a nested / regional air-quality domain driven
by lateral boundary conditions. It is the open-boundary counterpart of
`latlon3d_transport_cwc_3d_ppm`, which closes the same axes with a periodic zonal wrap and
polar caps.

## What it exercises

A mid-latitude patch — lon `[0, 60]`, lat `[20, 60]`, at **12×9×4** — so the zonal axis
does **not** close (real west/east walls) and the meridional axis does **not** reach the
poles (real south/north walls, face lengths `cos(lat) ~ O(1)`, not the vanishing polar
faces). Four rules fire at once:

- two bare-`D` **open-wall continuity** divergences, `face_flux_divergence_lon_open_bc` and
  `face_flux_divergence_lat_open_bc`, closing the air-mass budget with every wall face
  *differenced* (not wrapped, not omitted); and
- two compound **monotone PPM flux-form inflow** rules, `ppm_flux_D_lon_mono_inflow_bc` and
  `ppm_flux_D_lat_mono_inflow_bc`, carrying the tracer with a prescribed boundary halo.

The vertical axis is passive (no `lev` transport): this gate isolates the two new lon/lat
wall closures. The model carries the **air mass** `m` and the **tracer mass** `mq`, derives
`q = mq/m`, and supplies a prescribed lateral-boundary concentration `qbc_{w,e,s,n} == 1`
in a one-cell halo along each wall. The steady air-mass flux is positive across the patch,
so the **west and south walls are inflow** and the **east and north walls are outflow**.

## The property under test

The wall flux is built as a full-order one-sided PPM reconstruction over the prescribed
halo, upwind-selected by the sign of the wall-normal wind:
`F_wall = max(w,0)·(halo reconstruction) + min(w,0)·(interior reconstruction)`. On inflow
the halo is donated into the domain; on outflow the interior parabola leaves. This is
strictly more general than the 1-D `cartesian_uniform_1d/rules/ppm_D_inflow_bc`, which
hard-codes one inflow wall and one outflow wall for a one-signed speed.

At `q == 1` with the halo `== 1`, every CW84 slope is a difference of equal values and is
therefore exactly `0.0` in IEEE; the reconstruction of the uniform halo is exactly uniform,
so each open wall flux collapses to `F = M` **bitwise** at every face including the walls.
The PPM tracer divergence then telescopes onto the independent open-wall continuity twin to
the last bit, `dmq/dt` equals `dm/dt`, and `dev = mq − m` stays exactly `0`.

**Measured** (julia, `Tsit5`, reltol 1e-10 / abstol 1e-12): `dev` is exactly `0.0` at t=0
and t=0.1, `Linf`, against a 1e-12 atol. The walls are genuinely open and the flow is
non-degenerate — the air mass moves `max|Δm| ≈ 0.253` (~25% of a cell) over the span as
mass enters and leaves through the walls — so `dev = 0.0` is a real free-stream measurement,
not a quiescent no-op. Build+solve ≈ 42 s.

## Two things to know

**1. The reference is the constant zero, so `0.0` proves nothing by itself.** Non-vacuity
is established by `problems/_falsifier_regional_inflow.esm` — the same problem with all four
`qbc_{w,e,s,n}` halos set to `1.000001` instead of `1`. With a mismatched inflow
concentration the open-wall PPM flux no longer collapses to the mass flux, the tracer stops
tracking the air mass at the inflow walls, and `dev` reaches **1.872164546719e-6** at
t=0.1 — six orders above the 1e-12 tolerance — so the case correctly **fails**. This is what
makes the gate meaningful: it proves the boundary condition is *live*, not a silently-dead
closure. The falsifier is regenerated from this problem plus the single halo change, so the
two cannot drift apart; re-run it if you change this problem.

**2. Only the t>0 assertion discriminates.** `dev` is a state integrated from `dev(0) = 0`,
so its t=0 value is `0.0` by construction even under perturbation — a cheap sanity check on
the initial condition, not evidence about the scheme.

## Why a CWC gate rather than a convergence gate

The one-sided boundary *reconstruction* order is already characterized on the 1-D grid by
`tests/conformance/convergence/advection_1d_inflow_ppm` (same CW84 chain, same halo-fill
closure). What is new here, and what this gate pins, is **conservation and free-stream
preservation on the open regional domain** — the load-bearing property for a lateral
boundary condition and the one that fails catastrophically (a uniform tracer growing
structure at the walls) if the tracer and continuity wall closures disagree by even one bit.

## Mesh and cost

12×9×4. NLON=12 and NLAT=9 give both horizontal axes a genuine mid-domain interior beyond
the six boundary slabs each wall closure needs (three per wall, the upwind-biased 7-cell
support). NLEV=4 is a small passive vertical dimension. The property is bitwise exactness
and is resolution-independent; build cost tracks stencil width, not cell count.
