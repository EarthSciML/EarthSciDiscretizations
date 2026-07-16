# The monotone-PPM 3-D free-stream gate

Consistency with continuity (CWC) for the **production** reconstruction: the monotone
CW84 PPM scheme, carried in flux form on all three axes of the `latlon3d` grid at
**12x7x9**. The donor-cell twin `latlon3d_transport_cwc_3d` gates the same joint
agreement (zonal wrap, polar face omission, area weights, wall closure) for the
1-cell donor scheme; this case gates it for the scheme we actually intend to run.

## What it asserts

The model carries the air mass `m` alongside the tracer mass `mq` and derives
`q = mq/m`, so six rules fire at once: three bare-`D` face-flux-divergence rules
closing continuity (lon periodic, lat polar, lev no-flux) and three compound monotone
PPM flux-form advection rules carrying the tracer. At `q == 1` every CW84 slope is a
difference of equal values and is therefore exactly `0.0` in IEEE; the edge value
collapses to `0.5*(1+1) + (1/6)*(0-0) = 1` exactly; the limiter leaves `aL = aR = 1`,
so the parabola is the constant 1 and its integrated donor-region average is exactly 1.
The face flux is then `F = max(w,0)*1 + min(w,0)*1 = M` bitwise, each compound
divergence collapses onto its bare-`D` continuity twin to the last bit, and `q` stays
exactly 1.

This is a **different argument** from the donor case's
(`F = 0.5*(M*(1+1) - |M|*(1-1)) = M`). The donor's proof does not carry over to PPM
and must not be cited here.

**Measured, not predicted** (2026-07-16, julia binding, `Tsit5`, reltol 1e-10 /
abstol 1e-12): `dev` is exactly `0.0` at t=0 and t=0.1, Linf, against a 1e-12 atol.
The claim above is now a measurement.

## Two things to know before touching this case

**1. The assertion reference is the constant zero, so `0.0` proves nothing by itself.**
A dead gate passes identically to a live one. Non-vacuity is established by
`problems/_falsifier_ppm3d.esm` — the same problem with `ic(mq)` multiplied by
1.000001, so `q = 1+1e-6` — which makes `dev` drift at rate `(q-1)*div(M)`. Measured:
`dev` reaches **1.082573779618968e-6** at t=0.1 at 12x7x9 (and 1.0372392631220615e-6
at 7x7x7) — six orders above the atol — and the case correctly **fails**. It is
regenerated from this
problem plus the perturbation, so the two cannot drift apart. If you change this
problem, re-run the falsifier.

**2. Only the t=0.1 assertion discriminates.** `dev` is a state integrated from
`dev(0) = 0`, so the t=0 assertion reads `0.0` by construction — it holds even in the
falsifier. It is a cheap sanity check on the initial condition, not evidence about the
scheme.

## Cost

~370 s wall, ~3.5 GB peak RSS. The donor twin runs in under 2 minutes; the difference
is the monotone PPM build, not the solve. Under esm-spec 9.7.3 template bodies are
inlined at load and the Expression AST has no `let`/sharing node, so the reconstruction
chain (slopes -> limited slopes -> edge values -> limiter -> parabola -> flux)
re-inlines each stage into the next. The affine access-kernel build upstream is what
made the case runnable at all; before it, runs at 12x7x9 and 7x7x7 were both killed
before producing a result.

Build cost is set by the number of structural boxes, which tracks the **stencil width**
rather than the grid: 12x7x9 (756 cells) builds no slower than 7x7x7 (343 cells).
The mesh is 12x7x9 because NLAT=7 is the floor for the monotone lat rule (three cap
rows at each pole plus a non-empty interior) and the property under test is bitwise
exactness, which is resolution-independent — but 12 lon cells and 9 levels exercise a
genuine mid-domain interior on those two axes rather than the all-boundary degenerate
case that N=7 produces.
