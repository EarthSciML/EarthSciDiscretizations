# Repository Layout

EarthSciDiscretizations holds two things:

1. A Julia package (`src/`, `test/`) that implements discretization operators
   and the discretization pipeline.
2. A catalog of authoritative discretization **rule files** (`discretizations/`)
   that declaratively describe how continuous PDE operators map onto discrete
   stencils. Rule files are validated against the EarthSciSerialization
   discretization schema (§7) and executed by the ESS rule engine.

## Top-level layout

```
.
├── src/                  Julia package source
├── test/                 Julia package tests
├── docs/                 Documenter.jl documentation sources
├── discretizations/      Discretization rule JSON files (catalog)
│   ├── finite_difference/
│   ├── finite_volume/
│   └── spectral/
├── .github/workflows/    CI workflows (Tests, Documentation, …)
├── Project.toml          Julia package manifest
└── README.md
```

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
