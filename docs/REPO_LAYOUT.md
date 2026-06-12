# Repository Layout

EarthSciDiscretizations is a **multi-language monorepo**. It holds:

1. A Julia package at the repo root (`src/`, `test/`, `Project.toml`) that
   implements discretization operators and the discretization pipeline.
2. Sibling package trees for the Python, Rust, and TypeScript bindings
   (`python/`, `rust/`, `typescript/`) that implement the cross-binding
   grid accessor runtime defined in [`GRIDS_API.md`](GRIDS_API.md).
3. A catalog of authoritative discretization **rule files**
   (`discretizations/`) that declaratively describe how continuous PDE
   operators map onto discrete stencils. Rule files are validated against
   the EarthSciSerialization discretization schema (§7) and executed by
   the ESS rule engine.

## Top-level layout

```
.
├── src/                  Julia package source
├── test/                 Julia package tests
├── Project.toml          Julia package manifest
├── python/               Python binding (earthsci-discretizations)
│   ├── pyproject.toml
│   ├── src/earthsci_discretizations/
│   └── tests/
├── rust/                 Rust binding (earthsci_grids crate)
│   ├── Cargo.toml
│   ├── src/
│   └── tests/
├── typescript/           TypeScript binding (@earthsci/grids)
│   ├── package.json
│   ├── tsconfig.json
│   ├── src/
│   └── tests/
├── docs/                 Documenter.jl docs + cross-binding specs (GRIDS_API.md, …)
├── discretizations/      Discretization rule JSON files (catalog)
│   ├── finite_difference/
│   ├── finite_volume/
│   ├── grids/            Grid-family schemas + fixtures (cartesian, lat_lon, mpas, vertical, …)
│   └── spectral/
├── tests/                Cross-binding conformance fixtures (grids, rules, discretization)
├── tools/                Doc plot / rule-matrix rendering scripts
├── .github/workflows/    CI workflows (Tests, Python, Rust, TypeScript, …)
└── README.md
```

## Per-binding CI

Each language subtree has its own GitHub Actions workflow, scoped by path
filters so a PR touching only `python/` doesn't rebuild the Rust crate and
vice-versa:

| Binding     | Subtree         | Workflow                           |
|-------------|-----------------|------------------------------------|
| Julia       | `src/`, `test/` | `.github/workflows/Tests.yml`      |
| Python      | `python/`       | `.github/workflows/Python.yml`     |
| Rust        | `rust/`         | `.github/workflows/Rust.yml`       |
| TypeScript  | `typescript/`   | `.github/workflows/TypeScript.yml` |

See [`GRIDS_API.md`](GRIDS_API.md) for the cross-binding API contract every
generator (Julia / Python / Rust / TypeScript) must conform to.

## `discretizations/` convention

Each scheme-family subdirectory contains one JSON file per named scheme.
File contents conform to the ESS §7 schema. The per-family split is a
starting convention; it may be refined (e.g., sub-splits by equation class
or grid topology) as content lands.

Adding a new rule:

1. Create `discretizations/<family>/<scheme-name>.json`
2. Validate against the ESS schema (once the ESS rule engine lands, this
   will be a preflight CI check)
3. Add a rule-application test under `test/`

## Not in this repo

- The rule engine itself (lives in EarthSciSerialization)
- Model composition and solver orchestration (lives in EarthSciModels)
- The underlying PDE specification schema (lives in EarthSciSerialization)
