# DUO→MPAS Voronoi-dual topology value-invention conformance

Cross-binding byte-identity harness for the **MPAS Voronoi-dual topology** of the
DUO icosahedral primal, expressed as a declarative value-invention FAQ (bead
`esd-heg.2` / D1b). Each MPAS dual cell IS a primal vertex; this pins, for dual
levels 1/2/3, the three dual-topology arrays that replace the topology half of the
imperative `_duo_voronoi_dual` (`src/grids/mpas.jl`, steps 3-4):

| array             | what it is                                          | how it is produced |
|-------------------|-----------------------------------------------------|--------------------|
| `n_edges_on_cell` | valence of each dual cell (5 pentagon / 6 hexagon)  | primal-vertex degree (D1a `vertex_face` group-by) |
| `cells_on_cell`   | neighbour primal vertices, in CCW angular order     | tangent-plane `atan2` sort + bridge-vertex `equijoin` intersection |
| `edges_on_cell`   | shared primal edge id per neighbour (same slot)     | `skolem_edge` + `rank` (D1a canonical edge numbering) |

The dual requires a subdivided primal, so the levels are **1/2/3**
(`n_cells` = 42 / 162 / 642); every level has exactly 12 pentagonal cells (the
icosahedral seeds), the rest hexagonal.

## Artifacts

- **`../../../../../discretizations/grids/duo/rules/voronoi_dual_topology.esm`** —
  the declarative value-invention FAQ document (RFC `semiring-faq-unified-ir`
  §7.3). Its FAQ equations are the angular-ordering KEY — the tangent-plane
  `atan2` bearing of each incident face centroid; the per-cell stable sort by that
  angle, the consecutive-pair bridge-vertex intersection (`cells_on_cell`), and the
  `skolem_edge` rank (`edges_on_cell`) are documented relational/geometry
  derivations the engine materializes in the bridge (as D1a's `cell_neighbors`).
- **`src/topology_faq.jl` → `voronoi_dual_topology_faq`** — the thin evaluation
  bridge: it routes the bridge-vertex intersection and edge-id lookup through
  `EarthSciSerialization.Relational.equijoin` / the D1a edge `rank`, and performs
  the one geometric step (the angular `atan2` sort) directly, since §5.5.1 rule 1
  forbids floats in relational keys. No shadow relational engine (AGENTS.md
  single-pathway; GRIDS_API §4.3).
- **`golden.json`** — this directory. The binding-neutral, byte-identical canonical
  serialization of the three arrays for dual levels 1/2/3.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge.

## Indexing (binding-neutral)

Julia is the reference binding (1-based ids, `0` = padding past the valence). The
golden is **0-based**, converted at the harness boundary — identical to
`tests/conformance/grids/duo/topology/regenerate_golden.jl`. `cells_on_cell`
entries are primal vertex ids; `edges_on_cell` entries are the D1a **canonical**
(sorted) primal edge ids (`primal_topology_faq(...).edges`), replacing the
imperative `_build_connectivity` `Dict`-insertion numbering. Each cell's ring is
serialized **ragged** — the active entries in angular order, length
`n_edges_on_cell` — so the padding-to-`max_edges` convention is not pinned (it is
reconstructible from `n_edges_on_cell`).

`rank` emits dense ids in each binding's native base (`julia=1`, `rust=0`,
`python=0`); the conformance comparison normalizes to canonical 0-based via
`canonical = reported − base`.

## Why it is byte-identical across Julia / Rust / Python

The DUO primal mesh is the same deterministic icosahedron in every binding. The
relational half (the bridge-vertex intersection, the edge-id rank) is a pure
function of a defined total order over integer tuples — never hash or insertion
order — guaranteed identical across the three bindings by the ESS determinism
contract (`CONFORMANCE_SPEC.md` §5.5). The one geometric step (the tangent-plane
`atan2` angular sort) rides the geometry pipeline's fixed evaluation order: the
angles are computed by the same operations in every binding, so the `sortperm`
produces the same CCW ring. So each binding's `voronoi_dual_topology_faq` output
serializes to these identical bytes.

## What the test asserts (`test/test_duo_voronoi_dual_topology_faq.jl`)

For each dual level 1/2/3:
1. **Match imperative output** — `n_edges_on_cell` and `cells_on_cell` (the
   angularly-ordered neighbour ring of primal vertex ids) equal the imperative
   `_duo_voronoi_dual` output exactly; `max_edges` matches. This is the
   byte-identity contract for the dual topology.
2. **Internal consistency** — each `edges_on_cell` slot addresses the canonical
   primal edge `{cell, neighbour}`; the dual neighbour relation is symmetric; the
   padding past each valence is zero.
3. **Byte-identity** — the canonical 0-based serialization equals `golden.json`.

## Regenerate

```bash
julia --project=. tests/conformance/grids/duo/voronoi_dual_topology/regenerate_golden.jl
```

Regenerate only when the dual-topology contract legitimately changes; an
unexpected diff means the value-invention output drifted.
