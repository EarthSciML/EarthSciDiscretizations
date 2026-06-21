# MPAS-SCVT spherical-topology leaf — conformance golden

Bead `esd-e5m.3` (D3) of the declarative MPAS-SCVT mesh generator (epic
`esd-e5m`).

`golden.json` pins the **binding-neutral, 0-based canonical serialization of the
deterministic integer connectivity** the spherical-topology leaf emits — the
spherical Delaunay `faces` and the dual Voronoi `n_edges_on_cell` /
`cells_on_cell` / `edges_on_cell` / `vertices_on_cell` — produced by
`scvt_voronoi_connectivity` (`src/grids/mpas_scvt_topology.jl`): the irreducible
convex-hull spherical Delaunay (canonical executor the s2bindings.rs S2B FFI,
bead `s2b-s7b`) composed with the landed declarative `voronoi_dual_topology_faq`
(bead `esd-heg.2`).

Two **non-degenerate** seeds (no four cospherically coplanar generators, so the
Float64 hull is exact and the connectivity is the unique, order-independent
topology every binding — including the exact-predicate S2B FFI — must reproduce):

| seed | generators | triangles (`2n−4`) | cells |
|------|-----------:|-------------------:|-------|
| `octahedron` | 6 | 8 | all 4-neighbour squares |
| `icosahedral_level1` | 42 | 80 | 12 pentagons + 30 hexagons |

Indexing is 0-based binding-neutral; `faces` / `cells_on_cell` are generator
(cell) ids, `vertices_on_cell` are dual Voronoi vertex (Delaunay triangle) ids,
`edges_on_cell` are the canonical (sorted) primal-edge ids; rings are ragged to
`n_edges_on_cell`. Every binding's leaf output MUST serialize to these identical
bytes — the determinism contract (`CONFORMANCE_SPEC.md` §5.5 / §5.7;
[`../../../../../discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md`](../../../../../discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md)).

The circumcentre **geometry** is not pinned here — it rides the **tolerance**
contract (its unit-norm-times-`R` invariant is asserted to tolerance in
`test/test_mpas_scvt_topology_leaf.jl`), not the byte-identity determinism
contract that governs the integer connectivity.

The icosahedral-L1 `cells_on_cell` is additionally byte-identical to the
imperative `_duo_voronoi_dual` (the ρ ≡ 1 SCVT regression), asserted in the same
test.

Regenerate after an intentional, reviewed change to the leaf:

```
julia --project=. tests/conformance/grids/mpas/scvt/topology_leaf/regenerate_golden.jl
```
