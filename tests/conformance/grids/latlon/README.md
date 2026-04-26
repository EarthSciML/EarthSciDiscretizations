# Lat-lon cross-language conformance harness

Verifies that the grid-accessor bindings (Python, Rust, TypeScript — and
Julia once `dsc-1ts` lands) produce numerically equivalent accessor output
for the `lat_lon` family on a shared set of query points.

Scope follows the 2026-04-20 mayor correction (bead `dsc-suu`): the `.esm`
grid artifact is a **declarative config**, not a serialized geometry blob.
Cross-binding conformance therefore compares **accessor outputs**
(`cell_center`, `neighbors`, `metric_eval`, `cell_area`) at pinned query
points rather than serialized geometry bytes.

## Layout

- `fixtures.json` — fixture opts and query points, shared across bindings
- `regenerate_golden.py` — reference-golden regenerator
- `golden/<fixture>.json` — reference accessor outputs and rule-eval coefficients
- `../../../../python/tests/test_lat_lon_conformance.py` — Python accessor test
- `../../../../rust/tests/lat_lon_conformance.rs` — Rust accessor test
- `../../../../rust/tests/lat_lon_rule_conformance.rs` — Rust rule-eval test
- `../../../../typescript/tests/lat_lon.conformance.test.ts` — TypeScript test

## Reference binding

`docs/GRIDS_API.md` §4.3 pins **Julia as the reference binding** for ULP
ties. The Julia `lat_lon` runtime is not yet on `main` (tracked by bead
`dsc-1ts`), so the initial golden corpus is produced by the Python binding.

Downstream tests treat the golden as authoritative — they don't care which
binding emitted it. Once `dsc-1ts` lands:

1. Add a Julia conformance test under `test/test_lat_lon_conformance.jl`
   mirroring the three existing per-binding tests.
2. Regenerate the golden from Julia (or keep the Python script as a
   cross-check); the contract is unchanged.

## Fixtures

Three fixtures covering both variants:

| Name        | Variant            | Shape             | R (m)   | Notes                              |
|-------------|--------------------|-------------------|---------|------------------------------------|
| `small`     | `regular`          | nlon=16, nlat=8   | 1.0     | Unit sphere. 14 query points.      |
| `realistic` | `regular`          | nlon=360, nlat=180| 6.371e6 | Canonical 1° grid. 17 query points.|
| `reduced`   | `reduced_gaussian` | nlat=8, variable  | 1.0     | Per-row `nlon_per_row`. 20 query points covering every row width. |

Query points target row-boundary cases (first row, last row — where N/S
neighbors are `None` under `pole_policy="none"`), longitude wrap (i=0 and
i=nlon-1 on the same row), and interior cells. The `reduced` fixture adds
cross-row neighbor mapping between rows of differing `nlon`.

## Comparison protocol

Per `docs/GRIDS_API.md` §4, every binding computes at every query point:

- `cell_center(j, i)` → `(lon, lat)` in radians
- `neighbors(j, i)` → `{W, E, S, N}` each mapping to `(j', i')` or `None`
  (pole under `pole_policy="none"`)
- `metric_eval(name, j, i)` for `name ∈ {J, g_lonlon, g_latlat, g_lonlat,
  ginv_lonlon, ginv_latlat, ginv_lonlat, area}`

Results are compared against `golden/<fixture>.json`:

- **Neighbor indices**: exact equality. `null` in the golden matches a
  pole sentinel (`None` / `null`) from the binding.
- **Float fields**: relative tolerance `1e-14` per §4.2. Lat-lon geometry
  is closed-form; the only transcendental hit is `sin` on cell-center
  latitudes, which is bit-stable across Linux/macOS libm.

## Rule-evaluator parity

Each golden additionally carries a `rule_evals` array under the same
`tolerance.relative` budget. Each entry pins the per-query-point bindings
(`R`, `cos_lat`, `dlon`, `dlat` for the FD `latlon` family) and the
evaluated `coeff` value for every stencil entry of the named rule:

- `rule` / `rule_path` — discretization name and rule JSON path
- `stencil_selectors` — copy of the rule's stencil selector list (axis,
  offset) for cross-checking
- `bindings_per_qp` — per-query-point binding map; bindings derived from
  each binding's grid accessor must match these to the same tolerance
- `stencil_coeffs` — evaluated `coeff` per stencil entry per query point;
  each binding's ESS-AST evaluator must reproduce these floats

Rules currently covered:

- `centered_2nd_uniform_latlon`
  (`discretizations/finite_difference/centered_2nd_uniform_latlon.json`)

The Rust evaluator (`earthsci_grids::eval_coeff`) is exercised against
this corpus by `rust/tests/lat_lon_rule_conformance.rs`.

## Indexing

All query points and accessor outputs are `0`-indexed for row `j` and
column `i`. The three implemented bindings use 0-based indexing natively;
the Julia binding, once added, will convert at the boundary (matching the
existing cubed-sphere pattern).

Neighbors for the `reduced_gaussian` variant use nearest-center column
mapping between rows of differing `nlon` — each binding must implement the
same mapping (`floor((i + 0.5) / nlon_from * nlon_to)`) to match.

## Regenerating the golden

```bash
python3 tests/conformance/grids/latlon/regenerate_golden.py
```

Commit the updated `golden/*.json` files. All bindings must re-verify.

## Running

### Python
```bash
cd python && PYTHONPATH=src pytest tests/test_lat_lon_conformance.py
```

### Rust
```bash
cd rust && cargo test --test lat_lon_conformance
cd rust && cargo test --test lat_lon_rule_conformance
```

### TypeScript
```bash
cd typescript && npm test -- lat_lon.conformance
```

All three exiting 0 ⇒ bindings conform on this corpus. Once the Julia
binding (bead `dsc-1ts`) lands, the Julia test will round out the matrix.
