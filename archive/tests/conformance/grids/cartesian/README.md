# Cartesian cross-language conformance harness

Verifies the four grid-generator bindings (Julia, Python, Rust, TypeScript)
produce numerically equivalent accessor output on a shared set of query
points for the `cartesian` family (1D / 2D / 3D, uniform + non-uniform).

Scope follows the 2026-04-20 mayor correction in bead `dsc-g0s`: the
`.esm` grid artifact is a **declarative config**, not a serialized geometry
blob. Cross-binding conformance therefore compares **accessor outputs**
(`cell_centers`, `cell_widths`, `cell_volume`, `neighbors`, `metric_eval`)
at pinned query points rather than serialized geometry bytes.

## Layout

- `fixtures.json` — fixture opts and query points, shared across bindings
- `golden/<fixture>.json` — reference accessor outputs (Julia-generated)
- `../../../../test/test_cartesian_conformance.jl` — Julia conformance test
- `../../../../python/tests/test_cartesian_conformance.py` — Python test
- `../../../../rust/tests/cartesian_conformance.rs` — Rust test
- `../../../../typescript/tests/cartesian.conformance.test.ts` — TS test

## Fixtures

Four fixtures exercise the 1D/2D/3D and uniform/non-uniform paths from
`docs/GRIDS_API.md` §2 (`nx/ny/nz + extent` uniform path; `edges`
non-uniform path):

| Name            | ndim | Shape           | Spacing     | Notes                          |
|-----------------|------|-----------------|-------------|--------------------------------|
| `small_1d`      | 1    | `nx=8`          | uniform     | Unit interval. 5 query points. |
| `small_2d`      | 2    | `nx=4, ny=3`    | uniform     | 8 query points (corners + interior). |
| `realistic_3d`  | 3    | `16 x 8 x 4`    | uniform     | Anisotropic extent. 8 qps.     |
| `nonuniform_2d` | 2    | `4 x 4` via `edges` | non-uniform | Both axes explicitly specified. |

Query points target corners (all-boundary), interior cells, and edge cells
so the `neighbors` accessor is exercised at every boundary configuration.

## Comparison protocol (per `docs/GRIDS_API.md` §4)

Each binding computes, at every query point:

- `cell_centers(...idx)` → per-axis center coordinates
- `cell_widths(...idx)` → per-axis widths
- `cell_volume(...idx)` → cell measure (length / area / volume)
- `neighbors(...idx)` → axis-aligned face neighbors (boundary sides omitted)
- `metric_eval(name, ...idx)` for every metric defined at the fixture's
  `ndim`: `volume`, `jacobian`, `g`, plus `dx`/`dy`/`dz` and
  `face_area_x`/`_y`/`_z` per axis that exists

Results are compared against `golden/<fixture>.json`:

- **Integer fields** (neighbor indices, axis, side, cell counts): exact
  equality required.
- **Float fields**: relative tolerance `1e-14` per §4.2 default. Cartesian
  accessors are pure multiplications and linear interpolation over Float64
  extents — no transcendentals — so the tight default holds.

### Neighbor serialization

The golden emits neighbors as a sorted list of
`{"axis": int, "side": int, "index": [int,...]}` entries. Order is
axis-major, low side (`-1`) before high side (`+1`). Axis and index values
are 0-based (binding-neutral). Boundary faces are omitted (no sentinel).

## Indexing

All query points and accessor outputs in the golden are **0-indexed**
(panel / axis / cell indices). Julia's `CartesianGrid` accessors are
1-indexed; the Julia runner converts at the boundary so the golden JSON
remains binding-neutral.

## Regenerating the golden

The committed golden is Julia-produced (per the reference-binding rule in
§4.3). To regenerate after a deliberate algorithm change:

```bash
julia --project=. tests/conformance/grids/cartesian/regenerate_golden.jl
```

Commit the updated `golden/*.json` files. All four bindings must then
re-verify.

## Running

### Julia
```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["--filter", "Cartesian cross-language conformance"])'
```

### Python
```bash
cd python && pytest tests/test_cartesian_conformance.py
```

### Rust
```bash
cd rust && cargo test --test cartesian_conformance
```

### TypeScript
```bash
cd typescript && npm test -- cartesian.conformance
```

All four exit 0 ⇒ bindings conform on this corpus.
