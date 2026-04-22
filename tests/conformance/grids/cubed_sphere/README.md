# Cubed-sphere cross-language conformance harness

Verifies the four grid-generator bindings (Julia, Python, Rust, TypeScript)
produce numerically equivalent accessor output on a shared set of query points.

Scope follows the 2026-04-20 mayor correction in bead `dsc-suu`: the `.esm`
grid artifact is a **declarative config**, not a serialized geometry blob.
Cross-binding conformance therefore compares **accessor outputs**
(`cell_centers`, `neighbors`, `metric_eval`, `area`) at pinned query points
rather than serialized geometry bytes.

## Layout

- `fixtures.json` — fixture opts and query points, shared across bindings
- `golden/<fixture>.json` — reference accessor outputs (Julia-generated)
- `../../../../test/test_cubed_sphere_conformance.jl` — Julia conformance test
- `../../../../python/tests/test_cubed_sphere_conformance.py` — Python test
- `../../../../rust/tests/cubed_sphere_conformance.rs` — Rust test
- `../../../../typescript/tests/cubed_sphere.conformance.test.ts` — TS test

## Fixtures

Two fixtures per `docs/GRIDS_API.md` §4:

| Name        | Nc | R (m)   | Purpose                                 |
|-------------|----|---------|-----------------------------------------|
| `small`     | 4  | 1.0     | Unit sphere. Fast. 18 query points.     |
| `realistic` | 16 | 6.371e6 | Earth radius. 20 query points.          |

Query points target interior cells, panel corners (0,0 and Nc-1,Nc-1), and
panel-edge cells (to exercise cross-panel connectivity).

## Comparison protocol (per `docs/GRIDS_API.md` §4)

Each binding computes, at every query point:

- `cell_centers(p, i, j)` → `(lon, lat)`
- `neighbors(p, i, j)` → `{W, E, S, N}` mapping to `(p', i', j')`
- `metric_eval(name, p, i, j)` for `name ∈ {J, g_xixi, g_etaeta, g_xieta,
  ginv_xixi, ginv_etaeta, ginv_xieta, area}`

Results are compared against `golden/<fixture>.json`:

- **Integer fields** (neighbor indices, cell counts): exact equality required.
- **Float fields**: relative tolerance `1e-14` per §4.2. The gnomonic metric
  crosses `D2^(3/2)` which is a documented float-library boundary; the primary
  SHA-identity check is relaxed to the ULP-fallback across bindings.

## Indexing

All query points and accessor outputs are `0`-indexed for panel, i, j.

Julia's `CubedSphereGrid` is internally 1-indexed; the Julia runner converts
at the boundary so the golden JSON remains binding-neutral.

## Regenerating the golden

The committed golden is Julia-produced (per the reference-binding rule in
§4.3). To regenerate after a deliberate algorithm change:

```bash
julia --project=. tests/conformance/grids/cubed_sphere/regenerate_golden.jl
```

Commit the updated `golden/*.json` files. All four bindings must then
re-verify.

## Running

### Julia
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
(Filters to `:conformance` tag with `Pkg.test(test_args=["--tag=conformance"])`.)

### Python
```bash
cd python && pytest tests/test_cubed_sphere_conformance.py
```

### Rust
```bash
cd rust && cargo test --test cubed_sphere_conformance
```

### TypeScript
```bash
cd typescript && npm test -- cubed_sphere.conformance
```

All four exit 0 ⇒ bindings conform on this corpus.
