# DUO subdivision as a value-invention FAQ pass (esd-heg.3)

`subdivide_refine.esm` expresses **one level of icosahedral refinement** as a
semiring-FAQ value-invention pass (RFC `semiring-faq-unified-ir` §7.3), the
declarative replacement for the imperative `_midpoint!` + body of
`_subdivide_icosahedron` in [`src/grids/duo.jl`](../../../../src/grids/duo.jl).

The chain (the §7.3 skolem / distinct / rank pattern applied to mesh refinement):

1. **`edge_set`** — an index-set-producing `aggregate` (`semiring: bool_and_or`,
   `distinct: true`) that mints one midpoint vertex per unique undirected edge of
   the face→vertex connectivity, keyed by a canonical `skolem` over the sorted
   endpoint vertex ids.
2. **`edges`** — exposed as a `kind: "derived"` index set via `from_faq`.
3. **`edge_dense_id`** — `rank` to dense ids (the array-backend index).
4. **`mid_coord`** — the geometry FAQ: the normalized average of the two
   endpoints (`mid_sum / sqrt(mid_sq)`), on the unit sphere.
5. **`child_vert`** — joins each face to its three edge-midpoints (`edge_of_face`)
   and emits four sub-faces in the order `(a,ab,ca),(b,bc,ab),(c,ca,bc),(ab,bc,ca)`,
   the same order as the imperative step.

The document is a **value-free structural spec** (the connectivity/coordinates
are parameter *shapes*; concrete data is supplied by the harness), mirroring
`tests/valid/aggregate/edge_enumeration_area_eff.esm` in EarthSciSerialization.
Index sizes are the canonical level-0→1 instance (Nc=20, Nv=12, Ne=30, 4·Nc
children); the pass shape is identical at any level. It validates against the ESS
`esm-schema.json`.

## What "byte-identical; matches the imperative step" means here

The value-invention engine numbers invented vertices/edges in the **§5.7
canonical (sorted) order** — *never* the imperative `push!`/`Dict` order, which is
not cross-binding stable. So the FAQ output is the imperative *mesh* (same vertex
coordinates, same triangles, same winding) under a **deterministic, cross-binding
canonical labeling**. `fixtures/canonical/edge_set_level0.json` pins those exact
bytes (the contract every binding must reproduce).

## Proof

[`test/test_duo_subdivision_faq.jl`](../../../../test/test_duo_subdivision_faq.jl)
drives the **landed ESS engine primitives** (`Relational.skolem_edge` / `distinct`
/ `rank` for the value invention; `eval_coeff` for the geometry FAQ — no shadow
evaluator, per `AGENTS.md`) and proves a single FAQ pass reproduces the imperative
`_subdivide_icosahedron` step exactly for levels 0→1 and 1→2: identical vertex
coordinate set, identical triangles **position-for-position with winding
preserved**, plus the determinism (sorted, permutation/winding-independent) and
canonical-byte contract.

---

# DUO primal-cell geometry as a FAQ (esd-heg.6)

`primal_geometry.esm` expresses the **per-cell primal geometry** of the icosahedral
mesh as a semiring-FAQ geometry pass (RFC `semiring-faq-unified-ir` §8.1), the
declarative replacement for the per-cell geometry loop in `build_duo_grid`
(and `_spherical_triangle_area`) in
[`src/grids/duo.jl`](../../../../src/grids/duo.jl). Inputs (both `CONST`) are
`face_vert` (the cell→corner-vertex connectivity from the topology/subdivision
passes) and `vert_coord` (unit-sphere vertex coordinates); `R` is the sphere radius.

All outputs are `sum_product` aggregates over the cell set:

1. **`cell_cart`** — the cell-center **centroid**: the normalized mean of the three
   corner unit vectors (`cent_sum` reduces over the 3 corners, `cent_nrm2` over
   space, `cent_unit = cent_sum / sqrt(cent_nrm2)`, scaled by `R`).
2. **`lon` / `lat`** — geographic coordinates of the unit centroid via `atan2` /
   `asin` (the ESS P0 inverse-trig leaves, bead `ess-9x1`).
3. **`area`** — the spherical-triangle **area** by **L'Huilier's theorem** on the
   three great-circle arcs (`acos` of the corner dot products), spherical-excess
   form `4·atan(sqrt(tan(s/2)·tan((s−a)/2)·tan((s−b)/2)·tan((s−c)/2)))` scaled by `R²`.

Like `subdivide_refine.esm`, the document is a **value-free structural spec** and
validates against the ESS `esm-schema.json`. Index sizes are the canonical level-0
instance (Nc=20, Nv=12); the pass shape is identical at any level.

## What "byte-identical; matches the imperative step" means here

The FAQ mirrors the imperative float ops exactly — squares are products (`x*x`, not
`^`); the centroid sum, the squared norm, and the side-arc dot products fold in
corner/space order; the L'Huilier tan-product is left-folded. The **only**
intentional divergence is the removal of the imperative `clamp`/`max` clipping
guards (the acceptance asks for deterministic formulas, **no clipping**): on a valid
icosahedral mesh every `acos`/`asin` argument and the L'Huilier radicand already lie
in range, so clipping never fires and the result matches `build_duo_grid`
**bit-for-bit (0 ULP)**. `fixtures/canonical/primal_geometry_level0.json` pins the
exact level-0 Float64 values (area, centroid, lon/lat) every binding must reproduce.

## Proof

[`test/test_duo_primal_geometry_faq.jl`](../../../../test/test_duo_primal_geometry_faq.jl)
drives the **landed ESS engine** (`eval_coeff` — the single-pathway passthrough, no
shadow evaluator) on the geometry expressions and proves they reproduce the
imperative `build_duo_grid` area/centroid/lon/lat to **0 ULP** across levels 0, 1, 2,
plus the physical invariant (areas sum to 4πR²), the schema-valid declarative
document, and the cross-binding canonical-byte contract (reproduced bit-for-bit by
both the imperative grid and the FAQ).

---

# Build-time level-fold: seed + compose refine × level (esd-heg.4)

`subdivide_fold.esm` expresses the **full level-N mesh** declaratively as the base
icosahedron **seed** folded through the one-level `subdivide_refine` pass **N
times**, where `N = grid options.level` is a build-time constant. This is the
declarative replacement for the per-level driver loop `for _ in 1:level` of
`_subdivide_icosahedron` ([`src/grids/duo.jl`](../../../../src/grids/duo.jl)) — the
last recursive kernel.

1. **The seed** (the *literal/trivial FAQ constant*): the 12 golden-ratio
   (`φ = (1+√5)/2`) icosahedron vertices, pre-normalization, plus the 20 base
   faces. `subdivide_fold.esm` declares the **normalize-to-sphere FAQ**
   `seed_coord = vert_raw / sqrt(Σ_d vert_raw_d²)` (squares as products, 3-arg `+`
   left-assoc), so the declarative seed is **bit-for-bit identical** to the
   imperative `_icosahedron_vertices`.
2. **The fold**: each pass's output vertices + faces become the next pass's input
   parameters (`vert_coord`, `face_vert`). After `level` passes the level-N mesh is
   complete. The single pass is already proven equivalent to one imperative step
   (above), so the fold reproduces the imperative full mesh by induction.

## No native repeat affordance — caller-side unroll

ESS has **no native build-time repeat / iterate / fold construct** over a
build-time-constant level: there is no `repeat`/`iterate`/`level`/`loop` key in
`esm-schema.json`, and no whole-model value-invention runner that ingests
`index_sets` + `variables` + `equations` and materializes the state arrays (the
engine exposes the per-primitive functions `skolem_edge` / `distinct` / `rank` and
the expression evaluator, which callers compose). So the `level`-times unroll is
performed **caller-side**, exactly as the landed single-pass proof drives one pass.
`subdivide_fold.esm` therefore declares the seed (the new declarative content) and
documents the fold composition; the missing native affordance is filed as
upstream **EarthSciSerialization bead `ess-vnk`** (a build-time
*repeat-pass-`level`-times* engine primitive would let the fold itself be declared,
not just the seed).

## Byte-identity contract

The fold's `§5.7` canonical labeling — base vertices keep their seed ids, each pass
appends midpoint vertices in `distinct`+`rank` order — is deterministic and
cross-binding stable (the imperative `push!`/`Dict` order is not).
`fixtures/canonical/mesh_level{0,1,2}.json` pin the exact canonical cell→vertex
connectivity (`faces_canonical`) and the `icos_level{0,1,2}.esm` counts
(`20·4^L` cells); the test regenerates and compares them byte-for-byte.

## Proof

[`test/test_duo_subdivision_fold.jl`](../../../../test/test_duo_subdivision_fold.jl)
drives the same landed ESS primitives to build the declarative seed and fold the
refine pass `level` times, and proves the fold reproduces the imperative
`_subdivide_icosahedron` **full mesh** for levels 0, 1, 2: declarative seed
bit-identical to imperative, identical vertex coordinate set (coords to ULP),
identical triangles **position-for-position with winding preserved**,
`icos_level{0,1,2}` counts (`20·4^L` cells), the canonical full-mesh byte contract,
and seed-face permutation/winding-independence of the canonical labeling.

---

# DUO edge + MPAS Voronoi-dual geometry as a FAQ (esd-heg.7 / D2b)

`edge_dual_geometry.esm` expresses the **GEOMETRY half** of the imperative
`_duo_voronoi_dual` ([`src/grids/mpas.jl`](../../../../src/grids/mpas.jl), steps
1 / circumcenters / 5 / area) and the **DUO Tier-U** `edge_length` / `cell_distance`
([`src/grids/duo.jl`](../../../../src/grids/duo.jl)) as a semiring-FAQ geometry pass
(RFC `semiring-faq-unified-ir` §8.1). It is the **companion** of
`../rules/voronoi_dual_topology.esm` (D1b): that document produces the dual
**topology** and the angularly-ordered incident-face ring (`sorted_vertex_faces`);
this one **consumes** that ring (`ring_face` / `ring_face_next`) and produces every
dual geometric metric. The MPAS cell **is** the DUO vertex; the MPAS dual vertex
**is** the DUO face. Inputs (all `CONST`): `vert_coord` (recovered-unit cell
direction = `duo.vertices/R`), `centroid` (DUO face centroid = D2a `cell_cart/R`),
the `face_vert` / `edge_cell` / `edge_face` / `ring_face` connectivity, and `R`.

Outputs:

1. **`circ_x/y/z`** — the DUO face **circumcenter** (MPAS dual vertex): the
   normalized cross-product sum `a×b + b×c + c×a`. The imperative outward-orientation
   flip (`cc·centroid < 0 → negate`, which fires for **every** face on the
   consistently-wound mesh) is reproduced as the **standard cross sum × −1.0** — so
   even a zero component keeps the imperative's `−0.0` (reversing the operand order
   would give `x−x = +0.0` and break byte-identity on the symmetric faces).
2. **`dc_edge`** — cell-center-to-cell-center great-circle arc `R·acos(·)`; the same
   formula and inputs as the DUO Tier-U **`edge_length`** (vertex-to-vertex arc).
3. **`dv_edge`** — Voronoi vertex-to-vertex arc `R·acos(·)` over the two incident
   circumcenters.
4. **`cell_distance`** — the DUO Tier-U centroid-to-centroid arc (distinct from
   `dv_edge`, which uses circumcenters).
5. **`lon_cell` / `lat_cell`**, **`lon_edge` / `lat_edge`** — geographic coordinates
   of the cell centres and the normalized edge midpoints (`atan2` / `asin`).
6. **`area_cell`** — the Voronoi dual-cell area as the spherical-excess (L'Huilier)
   **fan** from the cell centre over consecutive circumcenter pairs in the D1b angular
   ring, summed and scaled by `R²`.

Like the other passes, the document is a **value-free structural spec** and
validates against the ESS `esm-schema.json`. Index sizes are the canonical level-0
dual instance (the icosahedron→dodecahedron dual: cells=12, faces=20, edges=30,
ring=5); the pass shape is identical at any level (the ring becomes ragged 5/6 once
subdivision introduces hexagons).

## What "byte-identical; matches the imperative step" means here

The FAQ mirrors the imperative float ops exactly — squares are products (`x*x`, not
`^`); dot products and the centroid/midpoint sums fold in space order; the L'Huilier
tan-product is left-folded; the circumcenter is the standard cross sum negated (the
universal flip, `−0.0` preserved). The **only** intentional divergence is dropping
the `clamp`/`max` clipping guards (the acceptance asks for deterministic formulas,
**no clipping**): on a valid icosahedral mesh every `acos` argument and the
L'Huilier radicand already lie in range, so clipping never fires and the result
matches `_duo_voronoi_dual` / `edge_length` / `cell_distance` **bit-for-bit**.
`fixtures/canonical/edge_dual_geometry_level1.json` pins the exact dual-level-1
Float64 values (circumcenter, `dc_edge` / `dv_edge` / `cell_distance`, `area_cell`,
cell/edge lon-lat) every binding must reproduce.

## Proof

[`test/test_duo_edge_dual_geometry_faq.jl`](../../../../test/test_duo_edge_dual_geometry_faq.jl)
drives the **landed ESS engine** (`eval_coeff` — the single-pathway passthrough, no
shadow evaluator) on the geometry expressions and proves they reproduce the
imperative `_duo_voronoi_dual` circumcenters / `dc_edge` / `dv_edge` / `area_cell` /
cell-and-edge lon-lat **plus** the DUO `edge_length` / `cell_distance` **bit-for-bit**
(strict Float64 bit identity, including `±0.0`) across dual levels 1, 2, 3 — consuming
the D1b `voronoi_dual_topology_faq(...).sorted_vertex_faces` ring for the area fan —
plus the physical invariant (dual areas sum to 4πR²), the schema-valid declarative
document, and the cross-binding canonical-byte contract (reproduced bit-for-bit by
both the imperative builder and the FAQ).
