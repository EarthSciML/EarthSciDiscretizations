---
title: "Implementing the intersection_area kernel across ESS bindings"
description: "Research note and implementation recommendation for the intersection_area geometric kernel-factor (RFC §8.1) across the Julia, Rust, and Python implementations of earthsci-toolkit. Surveys planar vs spherical polygon-clipping libraries and shared C++ engines, and resolves the cross-binding conformance problem that floating-point geometry makes unavoidable."
---

> **Status:** Research note / implementation recommendation. Companion to
> [*A semiring-parameterized FAQ IR for ESS arrayops*](semiring-faq-unified-ir.md)
> (§8.1, the required `intersection_area` op) and to
> [*Implementing the build-time relational engine across ESS bindings*](relational-engine-implementation.md)
> (§8, the conservative-regridding case study). **Target repo:** EarthSciSerialization.

---

## 1. Scope and the constraint that flips the calculus

`intersection_area(poly_a, poly_b)` returns the area of the geometric intersection of
two grid cells (RFC §8.1). It is a **kernel-factor leaf** — an opaque geometric
primitive the evaluator provides, not a semiring expression — used at **setup time**
to build first-order conservative regridding weights (`A_ij`, the
relational-engine note §8). Earth-science grids are on the **sphere** (lat-lon,
cubed-sphere), so the kernel must in general do **spherical** polygon clipping
(great-circle / parallel-meridian edges); treating lat-lon as a flat plane is wrong
near the poles and the antimeridian.

**The decisive difference from the relational engine.** The relational engine's
cross-binding determinism was *solvable bit-for-bit* by sorting integer tuples. This
kernel is the opposite: it produces a **floating-point area from polygon clipping**,
and FP clipping is irreducibly implementation-dependent — different libraries (or a
C++ engine vs. a pure-language port) differ in intersection ordering, area
summation order, robust-predicate strategy, and snapping. So:

> **Bit-for-bit identity across *independent* geometry implementations is
> unachievable. Conformance for `intersection_area` must be tolerance-based.**

This single fact inverts the shared-engine calculus. For the integer relational
engine a shared engine (DuckDB) was rejected because parallel-native + sort already
gave bit-identity for free. Here, parallel-native *cannot* give bit-identity, so a
**shared C++ core becomes genuinely attractive** — it collapses the problem from
"three algorithms" to "one algorithm, three builds," the only route to near-exact
agreement.

## 2. The architectural move that mostly dissolves the problem

`intersection_area` lives in the **static partition** (RFC §6.1), and its output is a
**serializable sparse weight matrix** (`ConservativeRegridding.jl`'s `Regridder` is
exactly this — relational-engine note §8). That means the cross-binding question can
be sidestepped:

> **Compute regridding weights *once*, with a single reference geometry
> implementation, serialize the sparse `A_ij` (+ `dst_areas`) artifact, and have the
> other bindings *load and apply* it.** The apply is the integer-indexed
> `sum_product` FAQ — bit-identical across bindings by the relational-engine spec.
> No binding except the reference needs a geometry kernel at all.

This is just the "materialize + cache the static partition" decision applied to
geometry. It is the **primary recommendation**: the hard, non-portable floating-point
work is done in one place; conformance then applies to the *apply* (bit-identical)
plus a one-time tolerance check on the serialized weights. Three independent
bit-identical geometry kernels are needed **only** if every binding must build
weights from raw grids with no shared artifact — a requirement worth avoiding.

The natural reference implementation is **Julia + GeometryOps.jl**, because
`ConservativeRegridding.jl` already uses exactly that stack for native spherical
regridding (relational-engine note §8). ESS would emit/consume its serialized
weights rather than re-deriving them three times.

## 3. Per-binding libraries (when an independent build *is* required)

If a binding must build weights itself (e.g. dynamic/online meshes with no
precompute step), here is the best available tool per language. Rust is the weak
point.

### Julia — GeometryOps.jl (native spherical; the reference)
`GeometryOps.jl` (0.1.x line, JuliaGeo) is the only mature Julia library doing
**native, non-approximate spherical** polygon intersection + area: `Spherical()`
manifold, `ConvexConvexSutherlandHodgman` clipping, `area(Spherical(), …)` via
Girard's theorem; `Planar()` and `Geodesic()` (WGS84, via Proj) also available. Pure
Julia, no C++ binary. It is what `ConservativeRegridding.jl` calls internally
(confirmed from source). **Use as the reference geometry implementation.** LibGEOS.jl
(GEOS, planar-only) and Meshes.jl (planar/convex-oriented clipping) are not the right
tools for the spherical case.

### Python — spherely (S2) primary; shapely/GEOS + projection fallback
- **`spherely`** (0.1.1, April 2026) — vectorized NumPy bindings to **Google S2** via
  `s2geography`; true spherical: `spherely.area(spherely.intersection(a, b), radius=R)`.
  The closest spherical mirror of the shapely API. **Caveat: pre-1.0, explicitly
  "early stage."**
- **`shapely`** (2.1.2, GEOS 3.13.1) — mature but **planar only**; usable as a
  spherical approximation by reprojecting each cell pair to a **local equal-area
  projection** (`pyproj`) and clipping in the plane.
- Precedent: the dominant Python conservative-regridding stack is **xESMF/ESMF**
  (SCRIP algorithm), which does true spherical overlap in **3D Cartesian** on the
  sphere and emits a sparse weight matrix — i.e. it already follows the
  "compute-weights-once, apply-many" pattern of §2.

### Rust — the forcing function: no native spherical clipper exists
There is **no robust native-Rust spherical polygon-clipping + area** library in 2026:
- `geo` (0.33) / `i_overlay` (7.x) — robust **planar** Boolean ops; `geo` *does* ship
  spherical/geodesic **area** traits (`ChamberlainDuquetteArea`, `GeodesicArea`/Karney),
  but its **intersection is planar only**. So you can measure a spherical polygon you
  already have, not clip two on the sphere.
- `georust/geos` (`geos` 11 / `geos-sys` 2, GEOS 3.14) — **planar**, same C++ core as
  shapely/LibGEOS; spherical only via per-cell local projection.
- `s2` (yjh0502, pure-Rust port of Go S2) — **polygon boolean ops are unimplemented**
  (`loop`/`polygon`/`shapeindex` pending); stalled at 0.0.x.
- `sphersgeo` (STScI) — has `.intersection()` + spherical area, but its docs warn it
  is **explicitly non-rigorous** (degenerate/touching cases unhandled); 0.1.x.

This is decisive: the "three independent native libraries" model **cannot work for
Rust at all**. A Rust binding that must build weights independently has only two real
options — **FFI to a shared spherical C++ core (S2)** (§4) or **consume the serialized
weights** (§2). Rust is therefore the strongest argument *both* for the §2
artifact-sharing approach *and*, where that's impossible, for a vendored shared S2
core rather than per-language libraries.

## 4. Shared C++ engine options

| Engine | Julia | Python | Rust | Same core ×3? | Geometry | Verdict |
|---|---|---|---|---|---|---|
| **GEOS** (3.14) | LibGEOS.jl | shapely | `georust/geos` | **Yes** (one `libgeos` C API) | **Planar** | Only true 3-language same-core option; planar ⇒ needs per-cell local projection for the sphere; deterministic (OverlayNG, robust predicates) and near-bit-identical *if the same GEOS version/build is pinned across all three* |
| **S2geometry** | ✗ (JSoC proposal only) | spherely / `s2geometry` (SWIG) | pure-Rust port (≠ core) | **No** (2026) | **Spherical** (correct) | Right math, but not one shared core across the three; real binding only in Python (+ R) |
| **Clipper2** | ✗ (no binding) | `pyclipr` | `clipper2` (FFI) + port | No | Planar (int64, most deterministic) | No Julia binding; area not a direct output |
| **CGAL** | ✗ | CGALPY (heavy) | ✗ | No | Planar, exact | Impractical to bind across three |

**Reading of the table:** the only engine that is genuinely one-core-across-three
*off the shelf* is **GEOS**, and it is **planar**. The only correct-spherical engine,
**S2**, is *not* an off-the-shelf 3-language shared core in 2026 (Python-only real
binding via spherely/s2geography; Julia uses GeometryOps' native re-implementation of
S2-quality logic; Rust has only an incomplete pure port). So no single *off-the-shelf*
engine is both shared across all three and natively spherical.

But the constraint is "off the shelf," not "possible." **S2geometry can be made the
shared core by vendoring the C++ and writing thin FFI bindings** — Julia via
CxxWrap/`ccall`, Python via pybind11 (or just reuse spherely), Rust via FFI to
`s2geography`. This is the *only* path that is simultaneously (a) natively spherical,
(b) capable of bit-identity, and (c) viable in Rust (which has no native spherical
clipper at all, §3). It costs binding work up front but resolves both the geometry
and the determinism axes at once. So the real choice is: **§2 (compute once, share
the artifact — give up nothing)** or, if independent per-binding builds are required,
**a vendored S2 shared core** — *not* three independent native libraries, which fail
on Rust and force tolerance-based conformance regardless.

## 5. Spherical vs planar — the accuracy axis

This is a correctness question independent of determinism, and the established
earth-science regridders are unanimous: **do the clip on the sphere, never on a flat
lat-lon plane** (poles and the antimeridian break the plane). Two equivalent
paradigms are in production use — **line-integral via Green/divergence theorem**
(SCRIP/Jones 1999; TempestRemap, with Gauss–Green) and **spherical polygon clipping +
spherical-triangle (excess) summation** (YAC/CDO; ESMF effectively). They agree at the
core (a triangle's spherical excess *is* the contour integral of its boundary).
`ConservativeRegridding.jl`'s `GeometryOps` path is the second paradigm (spherical
Sutherland–Hodgman + Girard area). ESMF/TempestRemap/YAC operate in **3D Cartesian on
the unit sphere**, which is what removes the pole singularity and antimeridian seam.

**The dominant error source is the edge model, and it is a real design decision here.**
A lat-lon cell edge running along a parallel is a *small circle*, **not a great
circle** — but S2, GeometryOps, ESMF (`GREAT_CIRCLE`) and TempestRemap all treat every
edge as a geodesic. Per the GMD 2024 "Truly conserving…" analysis, assuming
great-circle edges for a 30° lat-lon cell adjacent to the pole gives a **~4% area
error** (≈17% at 60° width, ≈1% at 15°); the fractional error scales with the **square
of the cell's longitude width**, so it is severe only for coarse polar cells and
**2+ orders of magnitude smaller at typical few-degree climate resolution**. Only
**YAC** natively distinguishes great-circle vs latitude-parallel edges per grid; the
standard mitigation elsewhere (XIOS) is to **densify** a parallel edge into many short
great-circle segments. So the `intersection_area` op's `manifold`/edge contract should
state the great-circle-edge assumption explicitly, and ESS should offer densification
for coarse lat-lon grids if polar accuracy matters.

The planar route — reproject each cell pair to a **local equal-area or gnomonic**
projection, clip in the plane (this is what a GEOS shared core would require) — is the
SCRIP-near-pole trick and is accurate only in the small-cell limit. For global meshes,
**prefer native spherical**; reserve planar-with-local-projection for the
GEOS-shared-core path or small-cell regional grids.

## 6. Conformance: tolerance-based, with a conservation invariant

Since exact cross-binding equality is unachievable (§1), the conformance suite for
`intersection_area` (and the weights built from it) must be **tolerance-based**:

1. **Combined relative + absolute tolerance per area:**
   `|a_x − a_ref| ≤ atol + rtol·a_ref`. Same shared core + pinned version/build →
   tight (`rtol ≈ 1e-9`). Independent spherical implementations (GeometryOps vs
   spherely) → an **empirically calibrated** `rtol`, plus a real `atol ≈ 1e-15·R²` to
   absorb **slivers**: near-tangent overlaps where the two clippers legitimately
   disagree on whether a tiny intersection even exists. Treat sub-`atol` areas as
   equal-to-zero ("present-but-tiny" and "absent" both pass) — that regime is where
   snapping/tie-breaking diverges and it does not affect weights.
2. **The physically meaningful gate is conservation, not per-cell agreement.**
   Make the primary conformance test the invariants — global mass conservation
   `Σ_j A_j·F_target[j] = Σ_i A_i·F_source[i]` and partition-of-unity `Σ_i W_ij = 1` —
   to a tight tolerance; per-pair `A_ij` equality is secondary, since it is the
   unstable sliver regime. Note a subtlety the precedent surfaces: first-order
   conservation is exact *only if computed cell areas equal true areas*, which edge
   approximations violate; established tools restore exact conservation with a
   **post-hoc global-mean/area correction** rather than perfect geometry, and the
   residual shrinks with resolution. `ConservativeRegridding.jl` sidesteps the
   normalization half of this by dividing by `dst_areas = Σ_i A_ij` (the row-sum of
   *computed* overlap areas), not the true target-cell area — so `Σ_i W_ij = 1` holds
   **by construction** regardless of edge error. ESS should follow that
   construction; conservation tolerance is then application-set and
   resolution-dependent, not a fixed epsilon.
3. **Declare the manifold.** The geometry interpretation (`Planar` / `Spherical` /
   `Geodesic`) is part of the op's contract (matching `ConservativeRegridding.jl`'s
   `Manifold`); two bindings can only be compared under the same manifold.
4. **GEOS as a cross-check oracle.** Even where it isn't the production kernel, a
   pinned GEOS (planar, local projection) is a useful third-party oracle for
   regression-testing the spherical kernels on small cells where planar ≈ spherical.

## 7. Recommendation

1. **Primary — single reference geometry impl + serialized weights (§2).** Compute
   `A_ij`/`dst_areas` once with **Julia + GeometryOps.jl** (native spherical; already
   the `ConservativeRegridding.jl` stack), serialize the sparse artifact, and have
   Rust/Python **load and apply** it. The apply is the bit-identical `sum_product`
   FAQ; no second/third geometry kernel needed. This makes Rust's lack of a spherical
   clipper a non-issue and reduces conformance to the apply + a one-time weight check.
2. **If a binding must build weights independently — vendor a shared S2 core, not
   three native libraries.** Rust has no native spherical clipper (§3), so the
   "three native libs" model is not even achievable; and independent FP clippers
   can't be bit-identical anyway. A **vendored S2geometry/s2geography C++ core with
   thin FFI bindings** (Julia CxxWrap, Python pybind11/spherely, Rust FFI) is the one
   route that is simultaneously spherical, bit-identity-capable, and Rust-viable. A
   lighter, tolerance-based fallback is acceptable on the Julia/Python side only —
   Julia → GeometryOps (`Spherical()`, the incumbent), Python → `spherely` — but
   **Rust still needs the shared/FFI core regardless**, which is why the shared S2
   core is the consistent choice once any independent build is required.
3. **Avoid the GEOS-planar route for global grids.** Pinning GEOS 3.14.x across all
   three (the only off-the-shelf shared core) gives bit-identity but only planar
   geometry + per-cell local projection — wrong for coarse/polar cells, and it
   diverges from the native-spherical regridder the Julia ecosystem already uses.
   Keep it only as a small-cell/regional fallback or a cross-check oracle.
4. **Op contract:** `intersection_area` carries a `manifold` flag; its conformance is
   declared tolerance-based (not bit-for-bit) in `CONFORMANCE_SPEC.md`, with the
   conservation/partition-of-unity invariants as the primary cross-binding tests.

The throughline: because the kernel is a static-partition, serializable factor, the
right design computes it **once** with the best spherical tool and shares the
artifact — turning an unsolvable three-way bit-identity problem into a solved
one-implementation problem, and leaving only the (bit-identical) apply to conform.

## 8. References

- **Verified source:** `JuliaGeo/ConservativeRegridding.jl` (GeometryOps-based,
  `Manifold`-selectable spherical/planar; STR-tree overlap; sparse `A_ij`),
  `JuliaGeo/GeometryOps.jl` (`Spherical()`/`Planar()`/`Geodesic()`, Girard area),
  `EarthSciML/EarthSciData.jl` (uses ConservativeRegridding for non-staggered grids).
- **Julia:** GeometryOps.jl <https://github.com/JuliaGeo/GeometryOps.jl>; LibGEOS.jl
  (GEOS 3.14, planar); Meshes.jl. No viable S2.jl (JSoC proposal only).
- **Python:** shapely 2.1.2 / GEOS 3.13.1 <https://github.com/shapely/shapely>;
  spherely 0.1.1 (S2 via s2geography) <https://github.com/benbovy/spherely>;
  s2geometry 0.14.0 (SWIG); xESMF/ESMF (SCRIP, 3D-Cartesian spherical, sparse weights)
  <https://xesmf.readthedocs.io/>.
- **Rust:** `geo` 0.33 (planar Boolean ops; ships geodesic *area* —
  `ChamberlainDuquetteArea`/`GeodesicArea` — but planar *intersection*) /
  `i_overlay` 7.x (planar engine behind `geo`); `georust/geos` 11 + `geos-sys` 2
  (GEOS, planar); `s2` (yjh0502, pure-Rust port — polygon boolean ops unimplemented);
  `sphersgeo` (STScI, spherical but explicitly non-rigorous intersection).
- **Shared engines:** GEOS 3.14.1 <https://github.com/libgeos/geos> (bindings:
  <https://libgeos.org/usage/bindings/>; OverlayNG determinism); S2geometry
  <https://github.com/google/s2geometry> + s2geography
  <https://github.com/paleolimbot/s2geography>; Clipper2
  <https://github.com/AngusJohnson/Clipper2>; CGAL 2D Boolean Set Operations.
- **Spherical-regridding precedent:** SCRIP/Jones 1999 (line integral, Lambert near
  poles); TempestRemap / Ullrich & Taylor 2015 (great-circle edges, Gauss–Green,
  overlap mesh) <https://github.com/ClimateGlobalChange/tempestremap>; ESMF
  (`GREAT_CIRCLE` vs `CART`, 3D-Cartesian) <https://earthsystemmodeling.org/regrid/>;
  YAC/Hanke et al. 2016 (search-and-clip, great-circle *and* latitude edges per grid);
  xESMF + sparselt (sparse-weight apply); "Truly conserving…" GMD 17:415 (2024) — edge
  model & 4–17% polar area error, post-hoc conservation correction; "Accurate and
  Robust Geometric Algorithms for Regridding on the Sphere" (EGUsphere 2026-636).
- Parent RFC: [`semiring-faq-unified-ir.md`](semiring-faq-unified-ir.md) §8.1;
  case study: [`relational-engine-implementation.md`](relational-engine-implementation.md) §8.
</content>
</invoke>
