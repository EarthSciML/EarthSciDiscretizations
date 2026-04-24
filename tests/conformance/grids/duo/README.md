# DUO cross-language conformance harness

Verifies the four DUO icosahedral-triangular grid bindings (Julia, Python,
Rust, TypeScript) produce numerically equivalent accessor output on a shared
set of query cells.

Scope follows the 2026-04-20 mayor correction on bead `dsc-gbo`: the `.esm`
grid artifact is a **declarative config**, not a serialized geometry blob.
Cross-binding conformance therefore compares **accessor outputs**
(`cell_centers`, `neighbors`, `metric_eval`) at pinned cell indices rather
than serialized geometry bytes.

## Layout

- `fixtures.json` — fixture opts and query cells, shared across bindings
- `golden/<fixture>.json` — reference accessor outputs (Julia-generated)
- `regenerate_golden.jl` — script that refreshes the golden from the Julia runtime
- `../../../../test/test_duo_conformance.jl` — Julia conformance test
- `../../../../python/tests/test_duo_conformance.py` — Python test
- `../../../../rust/tests/duo_conformance.rs` — Rust test
- `../../../../typescript/tests/duo.conformance.test.ts` — TS test

## Fixtures

Two fixtures per `docs/GRIDS_API.md` §4 (both `builtin://icosahedral/<level>`):

| Name        | level | R (m)   | Nc  | Purpose                                   |
|-------------|-------|---------|-----|-------------------------------------------|
| `small`     | 0     | 1.0     | 20  | Bare icosahedron on unit sphere (all 20). |
| `realistic` | 2     | 6.371e6 | 320 | Earth radius, 2× subdivided (20 cells).   |

`small` covers every cell. `realistic` samples the first 5 base-triangle
children of a few icosahedron faces — enough to exercise the midpoint cache
(cross-face connectivity) without bloating the golden.

## Comparison protocol (per `docs/GRIDS_API.md` §4)

Each binding computes, at every query cell `c`:

- `cell_centers(c)` → `(lon, lat)`
- `neighbors(c)` → `[nb0, nb1, nb2]` (0-based cell indices, `-1` = boundary)
- `metric_eval(name, c)` for `name ∈ {area, lon, lat, x, y, z}`

Results are compared against `golden/<fixture>.json`:

- **Integer fields** (neighbors, counts): exact equality required.
- **Float fields**: relative tolerance `1e-12` — a per-family override of
  §4.2's default `1e-14`. L'Huilier area accumulates acos/atan per triangle
  so libm variance across platforms lives in the `[1e-14, 1e-13]` band.
  Real algorithm drift lands well above `1e-12`.

## Indexing

All query cells and neighbor indices are 0-based. Julia's `DuoGrid` is
internally 1-indexed and uses `0` to mean "no neighbor"; the Julia runner
converts at the boundary so the golden JSON stays binding-neutral. Rust's
`neighbors(c)` returns `[Option<u32>; 3]`; the conformance harness maps
`None → -1`. Python / TypeScript already expose `-1` natively.

## Regenerating the golden

Julia is the reference binding per §4.3. To regenerate after a deliberate
algorithm change:

```bash
julia tests/conformance/grids/duo/regenerate_golden.jl
```

Commit the updated `golden/*.json` files. All four bindings must then
re-verify.

## Running

### Julia
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Python
```bash
cd python && pytest tests/test_duo_conformance.py
```

### Rust
```bash
cd rust && cargo test --test duo_conformance
```

### TypeScript
```bash
cd typescript && npm test -- duo.conformance
```

All four exit 0 ⇒ bindings conform on this corpus.
