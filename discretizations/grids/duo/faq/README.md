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
