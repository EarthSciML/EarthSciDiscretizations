# Regridding cross-binding conformance harness

Verifies the declarative regridder kernels (`regridding/<kernel>.esm`) satisfy
the RFC `pure-io-data-loaders` §8 regridding invariants **through the ESS
engine**, and pins the binding-neutral fixtures + golden.

Scope: bead `esd-47z.6` (the "ESD DRAIN" conformance bead). The per-kernel
**Julia evaluator** tests ship co-located with each `.esm`
(`regridding/{conservative_regrid_overlap_join,bspline_regrid,point_cell_average_regrid}_conformance_test.jl`,
beads `esd-9cb`/`esd-47z.3`/`.4`) and carry the rich end-to-end assertions; this
harness is the **shared cross-binding layer**.

## Layout

- `<kernel>/fixtures.json` — inputs (mesh / stations / sampled field), shared across bindings
- `<kernel>/golden.json` — reference outputs + invariant record
- `<kernel>/regenerate_golden.jl` — regenerates the golden
- `../../../test/test_regridding_conformance.jl` — Julia (reference binding) consumer

## Kernels & invariants (RFC §8)

| Kernel         | RFC §8 invariant                                                        | Cross-binding mode |
|----------------|-------------------------------------------------------------------------|--------------------|
| `bspline`      | reproduces the spline's polynomial degree exactly on staggered locations| **binding-independent** (closed-form oracle) |
| `point`        | cell mean of contributing stations; configured `missing_value` if empty | engine-bound (value-invention front-door) |
| `conservative` | `Σ_j A_j·F_tgt = Σ_i A_i·F_src` (conservation); `Σ_i W_ij = 1` (PoU)     | invariants binding-independent; per-pair areas tolerance-based (§5.8) |

## Cross-binding disposition (esd-47z.6)

The three kernels sit at three different points on the cross-binding feasibility
spectrum; this is the heart of the drain bead and is recorded here so the
boundary is explicit.

### `bspline` — binding-independent

A static-unroll polynomial weighted sum: `F_tgt[j] = Σ_k w_k(s[j])·F_src[base[j]+k]`,
no transcendentals, no geometry, no value-invention. The golden is the test
polynomial evaluated directly at the staggered targets (the **reproduction**
oracle). Every binding can reproduce it from the same `fixtures.json` (sample the
polynomial → build `F_src`/`base`/`s` → weighted sum). The Julia reference
consumer is **green**; Python/Rust/TS consumers are mechanically portable
**follow-up** ports. (`base` in the golden is 1-based, the ESS/Julia reference
index; a 0-based binding feeds `base − 1`.)

### `point` — engine-bound (value-invention)

The broad phase is a bin-Skolem equi-join materialised by the ESS **value-invention
front-door** (the `AbstractDict` evaluation path), which only the engine provides;
the bindings carry no shadow evaluator. So the **reference binding (Julia) is the
evaluating binding** for this kernel. What IS binding-neutral and shared: the
integer spatial-bin candidate set is byte-identical across bindings per
`CONFORMANCE_SPEC.md` §5.5 (no floats in keys), and the `fixtures.json` +
`golden.json` (analytic cell means + `missing_value`) are the shared contract a
future engine-backed binding consumes. The Julia consumer is **green**; a
cross-binding port is gated on a binding-level value-invention front-door (out of
scope for this bead — tracked as follow-up).

### `conservative` — invariants shared, geometry tolerance-based (ESS §5.8)

Per `CONFORMANCE_SPEC.md` §5.8 (the Conservative-Regridding Geometry Tolerance
Contract), the per-pair overlap areas `A_ij` are floating-point polygon-clip
results — GeometryOps (Julia) vs S2 (Python, Rust) — and are **tolerance-based,
deliberately NOT byte-identical**; ESS owns the per-binding clip kernels
(`ess-my4.4.*`). The **primary** cross-binding gate is the two physical
invariants, which are **exact by construction** and therefore binding-independent:

- **conservation** `Σ_j A_j·F_tgt = Σ_i A_i·F_src` — both sides from the same
  engine-computed areas, so exact regardless of edge-model error;
- **partition-of-unity** `Σ_i W_ij = 1` — normalised by `dst_areas = Σ_i A_ij`.

This harness's Julia consumer builds `A_ij` via the spherical clip, drives the
`.esm` assembly through the ESS engine, and asserts both invariants plus a
reference-binding `F_tgt`/`A_j` regression to `area_relative` (§5.8 tolerance
band). The cross-binding **geometry** conformance proper lives in ESS
(`tests/conformance/geometry/`, gate `ess-my4.4.8`); ESD reuses that kernel
(bead text: "conservative … existing kernel") rather than re-porting it.

## Regenerating

```
julia tests/conformance/regridding/bspline/regenerate_golden.jl        # closed form (JSON only)
julia tests/conformance/regridding/point/regenerate_golden.jl          # analytic (JSON only)
julia tests/conformance/regridding/conservative/regenerate_golden.jl   # needs ESS + GeometryOps
```

`bspline` and `point` regenerate from a throwaway `JSON`-only environment;
`conservative` activates a throwaway environment with `EarthSciSerialization` +
`GeometryOps` + `GeoInterface` to run the spherical clip.
