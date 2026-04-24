# Vertical cross-language conformance harness

Verifies that every binding (Python, TypeScript, Julia, Rust) produces
byte-identical canonical `.esm` bytes for the same vertical-grid options,
per `docs/GRIDS_API.md` §3.5 and §4.1.

## Why byte-identity (not ULP tolerance)

Unlike the cubed-sphere harness, every field in a vertical `.esm` is either:

- a pure rational (uniform sigma: `1 - k/nz`),
- a value the caller passed in verbatim (`z`, `theta`, `z_star` `levels`,
  explicit `eta` `ak`/`bk`), or
- a single IEEE-754 multiply-add (`eta` synthesized sigma `= ak/p0 + bk`).

There is no transcendental math, so the `docs/GRIDS_API.md` §4.2 ULP
fallback should not fire. The harness enforces `relative = 0.0` (strict byte
equality). If that ever breaks, investigate bindings before relaxing.

## Layout

- `fixtures.json` — per-fixture opts, tolerance, and the golden directory
  pointer.
- `../../../../discretizations/grids/vertical/*.esm` — the committed canonical
  `.esm` files. These serve as both the ESD fixture corpus **and** the
  conformance golden; storing a second copy under `golden/` here would leave
  two sources of truth for the same declarative config.
- `../../../../python/tests/test_vertical_conformance.py` — Python runner.
- `../../../../typescript/tests/vertical.conformance.test.ts` — TypeScript
  runner.
- `../../../../test/test_vertical_conformance.jl` — Julia runner.
- Rust runner lands with `dsc-31p`.

## Fixtures

Five fixtures covering all coordinate kinds specified in
`docs/GRIDS_API.md` §2.4, §3.2:

| Name                  | Coordinate | `nz` | Role      | Purpose                                      |
|-----------------------|------------|------|-----------|----------------------------------------------|
| `sigma_uniform_n16`   | `sigma`    | 16   | small     | Smallest resolution; pure rational levels.   |
| `sigma_uniform_n64`   | `sigma`    | 64   | realistic | Research-model-sized uniform sigma column.   |
| `z_troposphere_l32`   | `z`        | 32   | realistic | Geometric altitude, 0 m → 34 km.             |
| `eta_hybrid_l12`      | `eta`      | 12   | realistic | CAM-style hybrid sigma-pressure; ak/bk path. |
| `theta_isentropic_l10`| `theta`    | 10   | realistic | Isentropic surfaces 280 K → 380 K.           |

`sigma_uniform_n16` is the `small` fixture (per bead `dsc-bbj`). The
remaining four collectively constitute the `realistic` set: each is tied to
a distinct coordinate kind so the harness covers every code path in §2.4
without repeating the same math.

## Comparison protocol

For each fixture:

1. Binding calls `vertical(**opts)` / `grids.vertical(opts)` / … with the
   `opts` in `fixtures.json`.
2. Binding calls `to_esm()` / `toESM()` / … to obtain a declarative dict.
3. Both the binding output and the committed golden `.esm` have their
   `provenance.binding` and `provenance.binding_version` fields stripped
   (see `provenance_strip_fields` in `fixtures.json`). Everything else in
   `provenance` — `family`, `version`, `dtype`, `coordinate` — is required
   to match across bindings.
4. The stripped documents are serialized with `json.dumps(..., sort_keys=
   True, indent=2)` plus trailing newline (matching
   `discretizations/grids/vertical/regenerate_fixtures.py`) and compared
   byte-for-byte.

Byte equality is the whole test. There is no per-field ULP fallback.

## Provenance: why strip `binding` / `binding_version`

`provenance.binding` identifies the binding that emitted the `.esm`
(`"python"`, `"typescript"`, `"julia"`, `"rust"`) and `binding_version` is
its package version. Both are normative in the produced `.esm` (per §7) so
downstream tooling can tell who wrote the file. They are *not* supposed to
affect cross-binding equivalence, so the conformance canonicalization
erases them. All other provenance fields (family, version, coordinate,
dtype) must agree.

## Running

### Python

```bash
cd python && pytest tests/test_vertical_conformance.py
```

### TypeScript

```bash
cd typescript && npm test -- vertical.conformance
```

### Julia

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The runner lives at `test/test_vertical_conformance.jl` and carries the
`@testitem` tags `:grid`, `:vertical`, `:conformance` so the refinery can
filter it.

### Rust

Pending on `dsc-31p` (vertical accessor runtime for rust). The runner will
land at `rust/tests/vertical_conformance.rs` as part of that bead.

## Regenerating the golden

The "golden" IS `discretizations/grids/vertical/*.esm`. Do not regenerate
it from this harness. Regeneration happens via:

```bash
PYTHONPATH=python/src python3 discretizations/grids/vertical/regenerate_fixtures.py
```

After regenerating, every binding's conformance test must pass or the diff
is not landable. Per `docs/GRIDS_API.md` §4.3, Julia is the normative
reference binding for ULP ties; these fixtures have no ULP ties so Python
is the generator-of-record until the Julia binding lands and the committed
.esm bytes are re-emitted under its provenance (paper-trail only; the
stripped canonical form is unchanged).
