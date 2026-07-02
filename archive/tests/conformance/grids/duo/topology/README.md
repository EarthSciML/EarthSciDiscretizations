# DUO primal-topology value-invention conformance

Cross-binding byte-identity harness for the DUO **primal unstructured
connectivity** expressed as a declarative value-invention FAQ (bead `esd-heg.1`
/ D1a). It pins, for icosahedral levels 0/1/2, the four primal-topology arrays
the landed M3 relational engine produces from each level's `faces` connectivity:

| array            | what it is                                   | relational primitive |
|------------------|----------------------------------------------|----------------------|
| `edges`          | unique undirected vertex pairs (sorted)      | `skolem_edge` + `distinct` + `rank` |
| `edges_on_face`  | dense edge id of each triangle's 3 edges     | `rank`-id gather |
| `cell_neighbors` | cell across each edge of a face (`-1` = bdry)| `equijoin` (connectivity inversion) |
| `vertex_faces`   | sorted faces incident on each vertex         | `distinct` group-by |

## Artifacts

- **`../../../../../discretizations/grids/duo/rules/primal_topology.esm`** — the
  declarative value-invention FAQ document (RFC `semiring-faq-unified-ir` §7.3).
  Value-invents `edges` (the `edge_set` aggregate) and the `vertex_face_incidence`
  group-by basis; the two inversion joins (`edges_on_face`, `cell_neighbors`) are
  documented relational derivations the engine materializes.
- **`src/topology_faq.jl` → `primal_topology_faq`** — the thin evaluation bridge
  that routes edge enumeration / dedup / rank / inversion-join through
  `EarthSciSerialization.Relational`. ESD hosts no shadow relational engine
  (AGENTS.md single-pathway; GRIDS_API §4.3).
- **`golden.json`** — this directory. The binding-neutral, byte-identical
  canonical serialization of the four arrays for levels 0/1/2.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge.

## Indexing (binding-neutral)

Julia is the reference binding (1-based ids, `0` = no neighbor). The golden is
**0-based** with boundary sentinel **`-1`**, converted at the harness boundary —
identical to `tests/conformance/grids/duo/regenerate_golden.jl`. `edges` is the
ESS canonical undirected index set (`Relational.canonical_index_set_json`); the
dense arrays are compact JSON in the deterministic face / vertex order.

`rank` emits dense ids in each binding's native base (`julia=1`, `rust=0`,
`python=0`); the conformance comparison normalizes to canonical 0-based via
`canonical = reported − base`.

## Why it is byte-identical across Julia / Rust / Python

The `faces` connectivity is the same deterministic icosahedron in every binding,
and the value-invention itself (dedup order, dense numbering, skolem keys, join
order) is a pure function of a defined total order over tuples — never hash or
insertion order — guaranteed identical across the three bindings by the ESS
determinism contract (`CONFORMANCE_SPEC.md` §5.5; the relational + partition gate
landed in `ess-my4.3.11`). So each binding's `primal_topology_faq` output
serializes to these identical bytes.

## What the test asserts (`test/test_duo_topology_faq.jl`)

For each level 0/1/2:
1. **Match imperative output** — `cell_neighbors` and `vertex_faces` equal the
   imperative `_build_connectivity` / `_vertex_faces` output exactly (both are
   order-invariant); `edges` is the same undirected set, emitted in the canonical
   sorted order (replacing the imperative builders' binding-divergent
   `Dict`-insertion order).
2. **Internal consistency** — each `edges_on_face` slot addresses the correct
   sorted vertex pair, and neighbours across a slot share that edge id.
3. **Byte-identity** — the canonical 0-based serialization equals `golden.json`.

## Regenerate

```bash
julia --project=. tests/conformance/grids/duo/topology/regenerate_golden.jl
```

Regenerate only when the topology contract legitimately changes; an unexpected
diff means the value-invention output drifted.
