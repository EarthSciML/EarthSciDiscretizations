# Arakawa cross-language conformance harness

Verifies the four grid-generator bindings (Julia, Python, Rust, TypeScript)
produce numerically equivalent accessor output for the arakawa-staggered grid
family on a shared set of query points.

Scope follows the 2026-04-20 mayor correction (bead `dsc-suu`): the `.esm`
grid artifact is a **declarative config**, not a serialized geometry blob.
Cross-binding conformance therefore compares **accessor outputs**
(`cell_centers`, `u_face`, `v_face`, `corners`, `coord`, `neighbors`,
`metric_eval`) at pinned query points — plus the per-stagger
`variable_location` / `variable_shape` / `location_shape` tables, which are
integer-valued and must be identical bit-for-bit across bindings.

## Base-grid scope

The arakawa accessor runtime is parameterised over an `ArakawaBaseGrid`
abstraction. As of this harness, all four bindings implement `CartesianBase`;
an `ArakawaLatLonBase` wrapper that composes the lat-lon generator with an
arakawa transform has not landed in any binding (see
`discretizations/grids/arakawa/README.md` — the lat-lon fixtures committed
alongside the cartesian fixtures are declaration-only). **This harness
therefore exercises only the cartesian base.** When the lat-lon composition
lands, add a `lat_lon_*` fixture entry and extend `regenerate_golden.jl` to
build it.

## Layout

- `fixtures.json` — fixture opts, query points, tolerance. Shared across bindings.
- `golden/<fixture>.json` — reference accessor outputs (Julia-generated).
- `regenerate_golden.jl` — Julia script that (re)produces the golden files.
- `../../../../test/test_arakawa_conformance.jl` — Julia conformance test.
- `../../../../python/tests/test_arakawa_conformance.py` — Python test.
- `../../../../rust/tests/arakawa_conformance.rs` — Rust test.
- `../../../../typescript/tests/arakawa.conformance.test.ts` — TypeScript test.

## Fixtures

Two fixtures per `docs/GRIDS_API.md` §4 (size convention: small + realistic):

| Name        | Base      | nx×ny    | Extent                      | Purpose                              |
|-------------|-----------|----------|-----------------------------|--------------------------------------|
| `small`     | cartesian | 16×16    | `[0,1] × [0,1]`             | Unit square. Fast. Small query set.  |
| `realistic` | cartesian | 256×256  | `[-180,180] × [-90,90]`     | Geoscience-typical extent.           |

Query points per fixture target each location type (`cell_center`, `u_edge`,
`v_edge`, `corner`) at:
- domain corners (`(0, 0)`, `(nx, ny)` for corner/u_edge/v_edge ranges),
- interior points,
- asymmetric off-diagonal cells (to catch i/j swaps).

All five staggers (`A`, `B`, `C`, `D`, `E`) are exercised on every fixture —
though only their `variable_location` / `variable_shape` / `location_shape`
tables are stagger-dependent. The coordinate and neighbor tables are
stagger-independent (they depend only on the base grid), so they are sampled
once per fixture and checked against all bindings.

## Comparison protocol

At every query point, each binding computes:

- `cell_centers(i, j)` → `(x, y)` — checked against `coords.cell_center.xy`.
- `u_face(i, j)` under stagger C → `(x, y)` — checked against `coords.u_edge.xy`.
- `v_face(i, j)` under stagger C → `(x, y)` — checked against `coords.v_edge.xy`.
- `corners(i, j)` → `(x, y)` — checked against `coords.corner.xy`.
- `neighbors(loc, i, j)` → `{W, E, S, N}` each `(i, j)` or `null`/`None`/`nothing`.
- `metric_eval(name, i, j)` for `name ∈ {dx, dy, area}` at cell_center query points.

Per stagger, each binding also checks:

- `variable_location(var)` for `var ∈ {h, u, v}` — string equality.
- `variable_shape(var)` and `location_shape(loc)` — integer equality.
- `rotated` flag — boolean equality (true only for stagger E).

Comparison tolerances (per `docs/GRIDS_API.md` §4.2):

- **Integer fields** (neighbor indices, shapes): exact equality required.
- **Float fields**: relative tolerance `1e-14`. Cartesian base geometry is
  pure linear arithmetic (`xlo + (i + 0.5) * dx` and friends), so identical
  dispatch orders give bit-identical results across libm stacks; `1e-14` is
  a trivial margin that still flags real algorithm drift.

## Indexing

All query points and accessor outputs are `0`-indexed for `i`, `j`. Julia's
arakawa API is internally 1-indexed; the Julia runner converts at the
boundary so the golden JSON remains binding-neutral.

## How each binding routes location queries

Julia, Python, and TypeScript expose four named accessors
(`cell_centers`, `u_face`, `v_face`, `corners`) whose behaviour depends on
the stagger that labels the grid. To sample a specific location
independent of stagger, every test harness keeps two reference grids:

- `g_A` — stagger A (`u`, `v` both at cell_center). Used for `cell_centers`,
  `corners`, `neighbors(loc, ...)`, and `metric_eval`.
- `g_C` — stagger C (`u` at u_edge, `v` at v_edge). Used for `u_face`
  (routes to UEdge) and `v_face` (routes to VEdge).

Rust has both the named accessors and a generic `coord(loc, i, j)` helper;
the harness uses the named accessors for parity with the other three
bindings.

## Regenerating the golden

The committed golden is Julia-produced (per the reference-binding rule in
§4.3). To regenerate after a deliberate algorithm change:

```bash
julia --project=. tests/conformance/grids/arakawa/regenerate_golden.jl
```

Commit the updated `golden/*.json` files. All four bindings must then
re-verify.

## Running

### Julia
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
The conformance test is tagged `:conformance` and `:arakawa`; a
TestItemRunner filter can narrow the run to just those tags.

### Python
```bash
cd python && PYTHONPATH=src python3 -m pytest tests/test_arakawa_conformance.py
```

### Rust
```bash
cd rust && cargo test --test arakawa_conformance
```

### TypeScript
```bash
cd typescript && npm test -- arakawa.conformance
```

All four exit 0 ⇒ bindings conform on this corpus.
