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
