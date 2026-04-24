# MPAS cross-language conformance harness

Verifies the four grid-generator bindings (Julia, Python, Rust, TypeScript)
produce numerically equivalent accessor output on a shared set of query
points for the MPAS (unstructured Voronoi) grid family.

Scope follows the 2026-04-20 mayor correction (bead dsc-3nw): the `.esm`
grid artifact is a **declarative config**, not a serialized geometry blob.
Cross-binding conformance therefore compares **accessor outputs**
(`cell_centers`, `cell_center_cart`, `neighbors`, `cell_area`,
`edge_length`, `cell_distance`, `metric_eval`) at pinned query indices
rather than serialized geometry bytes.

## Layout

- `fixtures.json` — mesh arrays (inline), options, tolerance, and pinned
  query cells / edges, shared across bindings
- `golden/<fixture>.json` — reference accessor outputs (Python-generated,
  see below)
- `regenerate_golden.py` — Python regenerator; see "Regenerating" below
- `../../../../test/test_mpas_conformance.jl` — Julia conformance test
  (currently `@test_skip` pending bead dsc-7j0)
- `../../../../python/tests/test_mpas_conformance.py` — Python test
- `../../../../rust/tests/mpas_conformance.rs` — Rust test
- `../../../../typescript/tests/mpas.conformance.test.ts` — TypeScript test

## Mesh layout

MPAS meshes are unstructured and depend on a mesh file (typically NetCDF)
that is not bundled with this repo (see
`discretizations/grids/mpas/README.md`). The conformance harness therefore
carries a **small synthetic mesh inline** in `fixtures.json`: a
tetrahedral-Voronoi mesh with 4 triangular cells, 6 edges, and
`max_edges = 3`. Each cell has the other three as neighbors; the mesh
closes the sphere with total area `4 π R²`. This matches the
`_tetra_mesh` / `tetraMesh` fixture already used in
`python/tests/test_mpas.py` and `typescript/tests/mpas.test.ts`, keeping
the conformance corpus consistent with the per-binding unit tests.

## Fixtures

| Name              | Nc | R (m)     | Purpose                                           |
|-------------------|----|-----------|---------------------------------------------------|
| `tetra_small`     | 4  | 1.0       | Unit sphere. Fast. Exercises all 4 cells + 6 edges. |
| `tetra_realistic` | 4  | 6.371e6   | Earth radius. Same topology, real-scale magnitudes.  |

Both fixtures pin all 4 cells and all 6 edges as query indices, so every
accessor is exercised on every mesh element.

## Comparison protocol (per `docs/GRIDS_API.md` §4)

Each binding loads `fixtures.json`, constructs its native mesh-data
struct from the inline arrays, builds a grid via its `mpas` entry point,
and computes the following at every query cell:

- `cell_centers(c)` → `(lon, lat)`
- `cell_center_cart(c)` → `(x, y, z)`
- `neighbors(c)` → list of valid neighbor cell indices
- `cell_area(c)`
- `metric_eval(name, c)` for `name ∈ {lon, lat, area, x, y, z,
  n_edges_on_cell}`

And at every query edge:

- `edge_length(e)` (`= dv_edge`)
- `cell_distance(e)` (`= dc_edge`)
- `metric_eval(name, e)` for `name ∈ {lon_edge, lat_edge, dc_edge,
  dv_edge}`

Results are compared against `golden/<fixture>.json`:

- **Integer fields** (neighbor indices, `n_edges_on_cell`, dimensions):
  exact equality required.
- **Float fields**: relative tolerance `1e-12` per the schema's
  `tolerances.area / dc_edge / dv_edge` override. The mpas accessors do
  not pass through heavy transcendentals on the hot path (mesh arrays are
  stored, not recomputed), so empirical cross-binding drift on this
  corpus is well under `1e-14`; the tolerance mainly absorbs JSON float
  round-trip.

## Indexing

All cells and edges are `0`-indexed. The tetra mesh's adjacency arrays
use `-1` to mark "no neighbor" / "no edge" in the `cells_on_cell` /
`edges_on_cell` / `cells_on_edge` slots, matching the Python/Rust/TS
convention documented in each binding's `mpas` module.

## Reference binding (interim)

Per `docs/GRIDS_API.md` §4.3, Julia is the canonical reference binding
for ULP ties. The mpas-specific complication is that the Julia mpas
runtime is still landing on bead **dsc-7j0**
(`src/grids/mpas.jl`). Until that bead closes, **Python is the interim
reference binding** for the mpas family:

- `regenerate_golden.py` uses `earthsci_toolkit.grids.mpas` to emit
  `golden/*.json`.
- `test/test_mpas_conformance.jl` asserts the harness files exist and
  `@test_skip`s the accessor comparison with a clear message pointing at
  dsc-7j0. When the Julia `MpasGrid` lands, that skip guard flips and the
  Julia runner runs the same accessor comparison as the other three
  bindings.

The accessor output on this mesh is stored, not recomputed — so the
Python-vs-Julia reference choice does not materially affect observed
golden values beyond the existing `1e-12` tolerance. A follow-up polecat
can re-emit `golden/*.json` from Julia once dsc-7j0 lands; if the
numbers shift within tolerance, the harness will still pass for all
bindings without further changes.

## Regenerating the golden

```bash
python3 tests/conformance/grids/mpas/regenerate_golden.py
```

This reads the committed `fixtures.json` and rewrites
`golden/tetra_small.json` and `golden/tetra_realistic.json`. Commit the
updated files; all four binding tests must then re-verify.

To fully regenerate `fixtures.json` as well (only after a deliberate
mesh or query-set change):

```bash
python3 tests/conformance/grids/mpas/regenerate_golden.py --fixtures
```

The Python package must be importable: run with
`PYTHONPATH=python/src` from the repo root, or `pip install -e python/`.

## Running

### Julia
```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["--tag=mpas"])'
```
(Until dsc-7j0 lands, the accessor step is `@test_skip`; the structural
file-presence checks still run.)

### Python
```bash
cd python && PYTHONPATH="$(pwd)/src" pytest tests/test_mpas_conformance.py
```

### Rust
```bash
cd rust && cargo test --test mpas_conformance
```

### TypeScript
```bash
cd typescript && npm test -- mpas.conformance
```

All three currently-wired bindings exit 0 ⇒ bindings conform on this
corpus. The Julia runner will join once dsc-7j0 lands.
