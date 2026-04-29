# EarthSciDiscretizations — Grid API (Phase 0, normative)

Cross-binding public API contract for grid generators. This document is the
**normative** reference for Phases 1-7 of the ESD grid work. Every per-phase
implementation bead (28 total) must conform to the signatures, return-type
contract, and determinism guarantees defined here.

Status: **Phase 0, normative** — changes require a follow-up Phase 0 bead and
re-signoff from at least one polecat per binding.

Bindings in scope: **Julia, Python, Rust, TypeScript**. Go is intentionally
excluded as a generator (see §5); it remains a first-class *consumer* via the
`.esm` parser + validator contract in EarthSciSerialization.

Serialization contract: all grids round-trip through the EarthSciSerialization
(ESS) wire form, §5.4.6 canonicalization and §6 schema. This document
references those sections; it does not re-specify them.

---

## 0. Scope

**In scope:**

- Public function signatures for every binding's grid generators
- Return-type contract (in-memory grid object + `to_esm` / `toESM` lowering)
- Cross-language conformance protocol (byte-identical canonical `.esm`,
  float-ULP fallback)
- ESD's normative rule-evaluation policy: each binding's rule application
  and evaluation routes through the corresponding EarthSciSerialization
  (ESS) binding's official evaluator. ESD adds **no** numerical-evaluation
  code. Per-binding shadow evaluators are non-conforming. (§4.3)
- Ecosystem integration boundaries per binding
- Determinism requirements
- Per-family configuration schema shape (the fields; per-family parameter
  definitions land with their Phase N bead)

**Not in scope for this document:**

- Implementation of any specific grid family (Phases 1-7)
- Regridding / remapping (Phase 8)
- File-backed loaders such as MPAS `.nc` or DUO mesh files (Phases 6 / 7,
  though the loader *contract* is summarized in §10)
- Per-language project scaffolding — each binding's existing package layout
  is assumed to extend cleanly

---

## 1. Grid families

Every binding must expose the same set of named grid families. Families are
grouped into phases for implementation-bead scheduling, but the **API surface
is a flat namespace** in every binding.

| Phase | Family | Namespace leaf | Topology |
|---|---|---|---|
| 1 | Cartesian regular | `cartesian` | rectilinear box |
| 2 | Lat-lon (regular) | `lat_lon` | rectilinear sphere surface |
| 3 | Stretched lat-lon | `stretched_lat_lon` | rectilinear + pole/equator stretching |
| 4 | Cubed sphere (gnomonic) | `cubed_sphere` | 6-panel block-structured |
| 5 | Stretched cubed sphere (Schmidt / conformal) | `stretched_cubed_sphere` | 6-panel, non-uniform |
| 6 | MPAS (SCVT, loader-backed) | `mpas` | unstructured Voronoi |
| 7 | DUO / icosahedral triangulations | `duo` | unstructured triangular |

Phase 8 (regridding) and per-family variants (e.g., height-coordinate vertical
extensions) are out of scope for this document but will reuse §2's signature
conventions.

---

## 2. Public function signatures per binding

### 2.1 Naming contract

Every family named `<family>` in §1 is exposed as:

| Binding | Call form |
|---|---|
| Julia | `EarthSciDiscretizations.grids.<family>(; opts...)` |
| Python | `earthsci.grids.<family>(**opts)` |
| Rust | `earthsci_grids::<family>::builder().<opt>(…).build()` |
| TypeScript | `earthsci.grids.<family>(opts: <Family>Opts): Grid` |

Names **rhyme**: the function leaf, options keys, and option semantics are
identical across bindings. Idiomatic call forms differ; the *vocabulary* does
not.

### 2.2 Rules that apply to every binding

- **Option names are `snake_case`** in the wire form and in every binding's
  option key, *including* TypeScript and Rust. This is deliberate: cross-
  language diffs and log lines stay greppable. Rust's clippy `non_snake_case`
  allow is scoped to the generated options struct only.
- **No positional required parameters** in Julia / Python / TS. The only
  positional argument tolerated is a single options struct / dict / kwargs.
- **No global mutable state.** Every call is a pure function of its options.
- **No I/O on the hot path** unless the family is loader-backed (MPAS, DUO —
  see §10). Loader families take a `path` or `reader` option and the I/O
  contract is explicit.
- **Units are SI** unless a family explicitly declares otherwise in its
  schema (§8). Earth radius defaults to 6.371e6 m.
- **Default dtype**: `Float64` / `f64` / `number` / `np.float64`. Single
  precision is opt-in per binding; cross-language conformance (§4) is only
  promised at `Float64`.

### 2.3 Julia

```julia
using EarthSciDiscretizations
grid = EarthSciDiscretizations.grids.cubed_sphere(;
    Nc        = 48,
    R         = 6.371e6,
    dtype     = Float64,
    ghosts    = 3,
)
```

- Module path: `EarthSciDiscretizations.grids`, one function per family.
- All options are **keyword arguments**. Required options have no default; a
  missing required option raises `ArgumentError`.
- Returns a concrete struct `<Family>Grid{T} <: AbstractGrid` where `T` is
  the element type (`dtype`). See §3.1.

### 2.4 Python

```python
import earthsci
grid = earthsci.grids.cubed_sphere(
    Nc=48,
    R=6.371e6,
    dtype="float64",
    ghosts=3,
)
```

- Module path: `earthsci.grids`, one function per family.
- All options are **keyword-only** (declared via `*,` after the leading
  comma). Required options have no default; a missing required option
  raises `TypeError`.
- `dtype` is a string (`"float64"` | `"float32"`) or a `numpy.dtype`. String
  form is the canonical value for `.esm` round-trip.
- Returns an `earthsci.grids.Grid` — see §3.2.

### 2.5 Rust

```rust
use earthsci_grids::cubed_sphere;

let grid = cubed_sphere::builder()
    .nc(48)
    .r(6.371e6)
    .dtype(Dtype::F64)
    .ghosts(3)
    .build()?;
```

- Module path: crate `earthsci_grids`, submodule per family, `builder()` fn.
- Option setters on the builder use `snake_case` method names. A `.build()`
  call validates and returns `Result<Grid, GridError>`.
- Builders are `#[must_use]`; required options are enforced at `.build()`
  time (not at construction) so partial builders can be reused in tests.
- `dtype` is an enum `Dtype { F64, F32 }`. The default is `F64`.
- Returns a `Grid` value that implements the `earthsci_grids::Grid` trait.
  See §3.3.

### 2.6 TypeScript

```ts
import { earthsci } from "@earthsci/grids";

const grid = earthsci.grids.cubed_sphere({
  Nc: 48,
  R: 6.371e6,
  dtype: "float64",
  ghosts: 3,
});
```

- Module path: namespace `earthsci.grids`, one function per family.
- Single options-object argument typed as `<Family>Opts`. Required fields are
  non-optional in the TS type; missing required fields are compile errors.
- `dtype` is the string literal union `"float64" | "float32"`.
- Returns a `Grid` — see §3.4.

### 2.7 Option-name rhyme table (representative)

Every binding uses these keys (snake_case on the wire; Julia / Python / TS use
them verbatim; Rust converts to builder setters of the same name).

| Option | Type | Applies to | Meaning |
|---|---|---|---|
| `Nc` | int | cubed_sphere, stretched_cubed_sphere | cells per panel edge |
| `nx`, `ny`, `nz` | int | cartesian | cells per axis |
| `nlon`, `nlat` | int | lat_lon, stretched_lat_lon | cells in longitude / latitude |
| `R` | float | any spherical family | sphere radius (m) |
| `extent` | 2×3 float array | cartesian | axis-aligned bounding box |
| `ghosts` | int ≥ 0 | any | halo cell width |
| `dtype` | enum | any | element precision |
| `stretch` | struct | stretched_* families | per-family stretch config |
| `loader` | struct | mpas, duo | see §10 |

Per-family full schema lives in §8 + the family's Phase N bead.

---

## 3. Grid return-type contract

Every generator returns an in-memory **grid object**. Every grid object has a
lowering to a §6-schema-valid `.esm` value via a binding-specific method.
The lowered value, after `canonicalize()` per §5.4.6, is the cross-language
conformance ground truth (§4).

### 3.1 Julia

```julia
abstract type AbstractGrid end
# per family: struct CubedSphereGrid{T} <: AbstractGrid … end

to_esm(grid::AbstractGrid) :: Dict{String,Any}
```

- Grid structs wrap or interoperate with existing ecosystem types where
  natural (see §5). Structs **may** embed a `Meshes.Mesh` or a
  `ClimaCore.Spaces.AbstractSpace`, but the `to_esm` lowering is the
  canonical contract; ecosystem types are an implementation convenience.
- `to_esm(grid)` returns a plain `Dict{String,Any}` whose JSON lowering is
  §6-schema-valid.
- Julia is the comparison anchor for transcendental-ULP ties in
  **grid-object emission only** (§4.2). Julia is **not** a reference
  evaluator for rule semantics — every binding routes rule evaluation
  through its ESS counterpart (§4.3).

### 3.2 Python

```python
class Grid:
    # attributes documented per family
    def to_esm(self) -> dict: ...
```

- `earthsci.grids.Grid` is an `xarray.Dataset`-backed wrapper (subclass or
  composition — §5). The wrapper exposes grid-typed accessors (`.lon`,
  `.lat`, `.area`, …) and a `to_esm()` method returning a plain `dict`.
- `np.array` fallback: a pure-numpy grid object is available for environments
  without xarray, exposing the same `to_esm()` method.

### 3.3 Rust

```rust
pub trait Grid {
    fn to_esm(&self) -> serde_json::Value;
    fn family(&self) -> &'static str;
    fn dtype(&self) -> Dtype;
}
```

- Every family's returned struct implements `earthsci_grids::Grid`.
- `to_esm` returns a `serde_json::Value` whose string lowering (after
  canonicalize) matches the other bindings byte-for-byte.

### 3.4 TypeScript

```ts
interface Grid {
  readonly family: string;
  readonly dtype: "float64" | "float32";
  toESM(): object;
}
```

- Grid is a plain TS interface. `toESM()` returns a JS object suitable for
  `JSON.stringify` followed by §5.4.6 canonicalization.

### 3.5 Lowering invariant (normative)

For any family and any options `opts`:

```
canonicalize(to_esm_JL(grid_JL(opts)))
  == canonicalize(to_esm_PY(grid_PY(opts)))
  == canonicalize(to_esm_RS(grid_RS(opts)))
  == canonicalize(to_esm_TS(grid_TS(opts)))
```

where equality is byte-identical on the canonical wire form. The only
permitted deviation is the ULP-tolerance fallback defined in §4.

---

## 4. Cross-language conformance protocol

### 4.1 Primary check — byte-identical canonical form

Every binding emits §5.4.6-canonical `.esm` bytes. The cross-language
conformance test is:

```
sha256(canonicalize(to_esm_X(grid_X(opts))))
  == sha256(canonicalize(to_esm_Y(grid_Y(opts))))
```

for every pair `(X, Y)` of bindings and every opts vector in the conformance
corpus (maintained in `EarthSciSerialization/conformance/grids/`).

### 4.2 Float-ULP fallback (grid-object emission)

Some families (e.g., stretched cubed sphere) hit transcendental-function
boundaries where different math libraries disagree in the last few ULPs
of a *grid coordinate field*. When that happens, the primary SHA check on
grid emission fails and the fallback applies:

- Per-field relative tolerance: **1e-14** (default).
- Per-family override allowed; must be declared in §8 and reviewed by at
  least one polecat per binding.
- Tolerance is declared at the schema level, not per-test, so it is
  discoverable from `.esm` alone.
- For ULP ties on transcendental boundaries in **grid-object emission
  only**, Julia is the comparison anchor (libm-closest behavior on Linux
  and macOS; Julia's `@fastmath` discipline is off by default in the grid
  code). The conformance record carries the anchor binding name, the
  platform + math-lib fingerprint (per §6), and the max relative deviation
  observed per field.

This anchor applies **only** to byte/ULP comparison of grid-object output.
It does **not** authorize a "reference evaluator" for rule semantics. Rule
evaluation has no reference binding — see §4.3.

### 4.3 Cross-binding rule-evaluation regime — ESS evaluators only

ESD generates grids and authors rule libraries; it does **not** host
numerical-evaluation code. Each binding's rule application and evaluation
routes through the corresponding EarthSciSerialization (ESS) binding's
official simulation runner. This is normative across every ESD binding and
supersedes any prior "Julia is the reference evaluator" framing.

**1. Thin passthrough.** Each binding's `rule_eval` / `evaluator` /
`runner` entry point in ESD (where one exists at all — see §4.3.2) is a
thin passthrough to an ESS-binding official simulation runner:

| Binding | ESS official runners (representative) |
|---|---|
| Julia | `EarthSciSerialization.evaluate` (Symbolics / MTK), `tree_walk.jl` (large-system runtime) |
| Python | `earthsci_serialization.numpy_interpreter`, `earthsci_toolkit.simulate` |
| Rust | `earthsci_serialization::simulate`, `simulate_array` |
| TypeScript | the official ESS-TS evaluator entry point |
| Go | the official ESS-Go evaluator (consumer-only per §5.5) |

A binding **may** have more than one official ESS runner. Each must (a) be
documented as official in the binding's package docs, (b) consume the
ESS-canonical AST directly with no per-rule-shape dispatch in the runner
itself, and (c) be invokable as a public simulation API by users, not only
by tests. This list is illustrative; the binding's own AGENTS.md / README
is authoritative.

**2. No per-binding numerical-evaluation code in ESD.** ESD adds **no**
numerical interpreter, AST walker, stencil expander, "test-path
evaluator", "conformance regen" hand-walker, or shadow evaluator under
any binding. Rules are *authored* in ESD (as ESS-AST patches and rewrite
metadata); they are *applied* and *evaluated* by ESS. A PR that
introduces a separate numerical-evaluation pathway under any ESD binding
is non-conforming and must be rejected, regardless of how it is framed
("test only", "documentation only", "MMS convergence helper",
"validate against Julia", "regenerate goldens"). The anti-pattern list
in the workspace AGENTS / CLAUDE.md is canonical.

**3. Conformance fixtures live at the production-pipeline endpoint.** The
conformance fixture for any rule (or rule combination) is a pair of
ESS-canonical values:

- `input.esm` — the pre-discretization model state.
- `expected.esm` — the post-evaluation state.

The conformance test for binding *X* drives the canonical pipeline
end-to-end:

```
state_X    = ESS_X.load("input.esm")
applied_X  = ESS_X.discretize(state_X, esd_rules)   // ESD rules applied via ESS
result_X   = ESS_X.<official_runner>(applied_X)     // §4.3.1 runner
assert canonicalize(ESS_X.dump(result_X)) == bytes("expected.esm")
```

The "evaluator under test" is ESS's, not ESD's. ESD's role in the test is
the rule set passed to `discretize`, nothing more. There is no per-binding
stencil-walker or per-binding interpreter in the test path.

**4. Cross-binding agreement at the same endpoint.** Equivalence between
bindings is verified at `result_X` vs `result_Y` after ESS-canonical
serialization — not at any per-binding intermediate (stencil walker,
interpreter, or numerical helper). Julia and every other binding are
equally subordinate to ESS; no binding is a reference evaluator for any
other. The §4.2 ULP fallback applies only if the byte-identical check
fails on the same transcendental boundary that already affects grid
emission, bounded by the same per-field tolerance.

#### What this section replaces

This section retires the prior "Julia is the reference evaluator" regime
that previously lived here. That framing — originally a narrow ULP-tie
anchor for grid-object emission (now scoped to §4.2) — was being read as
authorization for per-binding shadow evaluators in ESD on the theory that
non-Julia bindings needed local numerical machinery to validate against
Julia output. They do not. ESS provides the evaluator in every binding;
ESD piggybacks. The regime above is the only one in force.

### 4.4 Conformance corpus

The conformance corpus is a tree of `(input.esm, expected.esm)` fixture
pairs maintained alongside each rule's definition (e.g.
`discretizations/<group>/<rule>/fixtures/<scenario>/{input,expected}.esm`).
Every implementation bead adds at least one fixture pair. Every binding's
CI runs the corpus through the §4.3 pipeline; failures block the merge
queue.

Fixtures are full ESS values, not `(family, opts) → SHA` tuples. The
SHA-only form pre-supposes a single canonical evaluator output, which is
exactly the regime §4.3 retires. Generation of `expected.esm` MUST come
from the canonical pipeline — driving `discretize` + an official ESS
runner — never from a hand-walked stencil or a side-channel evaluator.

---

## 5. Ecosystem integration per binding

Each binding has an ecosystem-integration hierarchy: **primary** types the
binding natively interoperates with, and a **fallback** for environments
without those dependencies. The primary layer is where grid structs embed
or return ecosystem objects; the fallback is always available so that the
binding has zero hard deps beyond its own package.

### 5.1 Julia

- **Primary**: `Meshes.jl` and `ClimaCore.jl`.
  - Cubed-sphere, lat-lon, and cartesian families expose a `Meshes.Mesh`
    view via `as_meshes(grid)`.
  - Cubed-sphere and lat-lon families expose a `ClimaCore.Spaces` view via
    `as_climacore(grid)` where applicable.
- **Fallback**: `StaticArrays` for geometry primitives. No heavy deps;
  every grid struct is usable without Meshes.jl or ClimaCore.jl installed.

### 5.2 Python

- **Primary**: `xarray` + `pyproj` + output compatible with `xESMF`'s grid
  input conventions. Grid wrappers are `xarray.Dataset`-backed.
- **Fallback**: pure `numpy.ndarray` grid object. Same `to_esm` lowering.
  `xarray`, `pyproj`, and `xESMF` are declared as optional extras.

### 5.3 Rust

- **Primary**: `geo` (geo-rust) for 2D geometric primitives; `proj`
  bindings for geographic projections when a family needs them.
- **Fallback**: `nalgebra` for linear-algebra primitives. `geo` and `proj`
  are optional features; `nalgebra` is the default.

### 5.4 TypeScript

- **Primary**: `d3-geo` for standard geographic projections.
- **Fallback**: hand-rolled spherical / curvilinear utilities inside
  `@earthsci/grids`. The TS ecosystem for climate-scale grids is thinner
  than the other three; we acknowledge this rather than force-fit an
  inappropriate dependency.

### 5.5 Go (exclusion — normative)

Go is **not** a grid-generator binding. Go has no first-class climate-grid
library, and building one would duplicate work better done in the other
four languages.

Go **is** a first-class **consumer**: the Go bindings in
EarthSciSerialization implement the §6 parser and §5.4.6 canonicalization
validator and can round-trip any `.esm` grid emitted by Julia / Python /
Rust / TypeScript. Go users consume grids produced elsewhere (e.g., via a
build-time CLI) and ship `.esm` files as data artifacts.

Any future proposal to add a Go generator must:

1. Point at a viable Go ecosystem (`go-geom`, `orb`, or equivalent with
   projection support).
2. Demonstrate the conformance corpus (§4.4) passes on Go.
3. File a new Phase 0 bead amending this document.

Until those three conditions are met, Go generator work is out of scope.

---

## 6. Determinism requirements

### 6.1 Within a (binding, platform, compiler, math-lib) tuple

Given the same options and the same build fingerprint, every call produces
**byte-identical** uncanonicalized output across runs. No PRNG-seeded
initialization, no insertion-order dependence in maps, no time or hostname
leakage into the payload.

### 6.2 Across platforms (same binding)

Byte-identical after `canonicalize()` (§5.4.6). Platform-dependent field
ordering in JSON is the main risk; canonicalization eliminates it.

### 6.3 Across bindings

Byte-identical after `canonicalize()` in the common case. Within
per-field ULP tolerance (§4.2) when transcendental-math divergence fires
on grid-object emission; Julia is the comparison anchor for those ULP
ties (§4.2). Rule evaluation across bindings has **no** reference
binding — every binding routes through its ESS official runner and
agreement is verified at the pipeline endpoint (§4.3).

### 6.4 Build fingerprint

Every `.esm` grid emitted by a generator carries a `provenance` block
(§6-schema-defined) with:

- Binding name + version
- Platform triple
- Compiler / runtime version
- Math library fingerprint (e.g., `libm-glibc-2.38`, `msvc-14.39`)
- Source SHA of the generator

This lets conformance failures be attributed to a specific build rather
than to "numerical flakiness."

---

## 7. Per-binding grid-object shape (normative fields)

Every grid object, regardless of binding, exposes the following fields.
Binding-specific accessors are named per §2.2 (snake_case). Phase N beads
add family-specific fields; the fields below are the **common minimum** that
every family provides.

| Field | Type | Shape | Meaning |
|---|---|---|---|
| `family` | string | scalar | e.g., `"cubed_sphere"` |
| `dtype` | string | scalar | `"float64"` \| `"float32"` |
| `topology` | string | scalar | `"rectilinear"` \| `"block_structured"` \| `"unstructured"` |
| `ghosts` | int | scalar | halo width (0 if unsupported) |
| `n_cells` | int | scalar | total cell count (interior + ghosts consistent with `ghosts`) |
| `lon`, `lat` | float | family-shaped | cell-center geographic coords (spherical families) |
| `x`, `y`, `z` | float | family-shaped | cell-center cartesian coords (cartesian family) |
| `area` | float | per-cell | cell face area (spherical) or volume (3D cartesian) |
| `provenance` | struct | scalar | §6.4 fingerprint |

Families add: staggered edges, panel connectivity (cubed_sphere), metric
tensors (curvilinear), adjacency lists (unstructured). Those land in the
family's §8 schema.

---

## 8. Configuration schema shape per grid family

Every family's config schema is a JSON-schema-subset document stored under
`discretizations/grids/<family>.schema.json`. Per-family content lands with
the family's Phase N bead. The **shape** is fixed here:

```jsonc
{
  "family": "<family>",
  "version": "1.0.0",
  "required": ["<...>"],
  "options": {
    "<option>": {
      "type": "int|float|string|enum|struct",
      "default": "<value or null if required>",
      "units": "m|rad|1|...",
      "doc": "<one-sentence description>"
    }
  },
  "tolerances": {
    "default_rel": 1e-14,
    "per_field": { "<field>": 1e-12 }
  },
  "ecosystem_hints": {
    "julia":  { "primary": "Meshes.jl",   "fallback": "StaticArrays" },
    "python": { "primary": "xarray",      "fallback": "numpy" },
    "rust":   { "primary": "geo-rust",    "fallback": "nalgebra" },
    "ts":     { "primary": "d3-geo",      "fallback": "hand-rolled" }
  }
}
```

The `tolerances.default_rel` field is §4.2's ULP-fallback knob; the
`per_field` overrides let specific families relax it only where needed.

Required common options (§2.7) that appear for a family must be declared
with their family-specific defaults; the **names** stay canonical.

---

## 9. Error contract

Each binding raises binding-idiomatic errors; the **conditions** under
which errors fire are shared.

| Condition | Julia | Python | Rust | TS |
|---|---|---|---|---|
| Missing required option | `ArgumentError` | `TypeError` | `GridError::MissingOption` | `TypeError` (compile-time preferred) |
| Invalid option value | `DomainError` | `ValueError` | `GridError::InvalidOption` | `RangeError` |
| Schema violation in output | `AssertionError` | `AssertionError` | `GridError::SchemaViolation` | `Error` |
| Loader-backed file not found | `SystemError` | `FileNotFoundError` | `GridError::Io` | `Error` |
| Conformance tolerance exceeded | `ErrorException` | `AssertionError` | `GridError::Conformance` | `Error` |

All errors carry a machine-readable `code` string (same across bindings),
so CI and conformance tooling can route on code rather than string match.

---

## 10. Loader contract (MPAS, DUO)

Loader-backed families take a `loader` option. Its shape is normative:

```
loader = {
  path:   <string>,     # filesystem path (may be sandboxed)
  reader: <enum>,       # "auto" | "nc4" | "mpas_mesh" | "duo_mesh"
  check:  <enum>,       # "strict" | "lenient"
}
```

- All I/O happens inside the generator call; the returned grid is fully
  materialized (no lazy readers crossing the API boundary).
- The generator must record the loader fingerprint (file SHA, reader
  version) into the grid's `provenance` block (§6.4).
- Loader errors raise binding-appropriate errors per §9.

Full file-format specs (`.nc` variant, `.mpas`, `.duo`) live in the
EarthSciSerialization repo and are referenced by hash here rather than
re-specified.

---

## 11. Versioning

This document is versioned with a `(major, minor, patch)` semver triple
embedded in the header of every `.esm` grid:

- **Major** bump: breaking change to signatures, return types, or
  conformance protocol. Triggers a new Phase 0 bead.
- **Minor** bump: additive (new family, new option with default, new
  ecosystem-integration primary). Existing corpus must still pass.
- **Patch** bump: documentation-only, ecosystem-hint refinements, or
  typo fixes.

Current version: **1.0.0** (first normative release, Phase 0).

---

## 12. Signoff

This document becomes binding once a polecat from **each** of the four
generator bindings (Julia, Python, Rust, TypeScript) has reviewed via
bead comments, per the Phase 0 acceptance criteria.

- [ ] Julia signoff
- [ ] Python signoff
- [ ] Rust signoff
- [ ] TypeScript signoff
- [ ] Go (consumer-only) acknowledged §5.5

Until all four signoffs land, downstream Phase 1 beads may draft against
this document but must not be marked landed.

---

## Cross-references

- EarthSciSerialization §5.4.6 — canonical wire form
- EarthSciSerialization §6 — grid schema
- EarthSciSerialization §7 — rule schema
- EarthSciSerialization per-binding AGENTS.md — list of official
  simulation runners per binding (authoritative for §4.3.1)
- ESD / workspace `CLAUDE.md` — Simulation Pathway rule and anti-pattern
  list (canonical for §4.3.2)
- ESD `docs/REPO_LAYOUT.md` — repo layout
- ESD `docs/rule-catalog.md` — Phase 0 rule inventory
