# EarthSciDiscretizations

Authoritative discretization rule files — and the conformance machinery that
keeps them honest — for the EarthSci model stack.

This repository has three roles:

1. **A catalog of declarative discretization rule JSON files** under
   [`discretizations/`](discretizations/) (finite-difference and
   finite-volume stencils, reconstructions, flux forms, and boundary
   conditions), validated against the EarthSciSerialization (ESS) §7
   discretization schema and executed by the ESS rule engine via the
   canonical `discretize → ArrayOp → eval` pipeline. The package side
   contributes structural lowering between rule forms
   (`src/stencil_lowering.jl`) and thin coefficient-evaluation
   passthroughs (`src/rule_eval.jl`). This catalog is the repository's
   primary, durable role.
2. **Grid runtimes and conformance machinery.** Grid constructors for
   cartesian, lat-lon, vertical, Arakawa-staggered, cubed-sphere, MPAS,
   and DUO topologies (`src/grids/`), mirrored by thin Python, Rust, and
   TypeScript bindings (`python/`, `rust/`, `typescript/`), plus the
   cross-binding conformance suite: golden fixtures under
   [`tests/conformance/`](tests/conformance/), the dynamic rule walker
   (`test/walk_esd_tests.jl`), and its machine-readable per-rule skip
   ledger (`test/junit-esd.xml`). One catalog, four bindings, one
   canonical pipeline — conformance-tested.
3. **Transitional reference operators.** Hand-coded finite-volume
   operators on cubed-sphere grids (`transport_2d`, `flux_1d`, PPM
   transport/reconstruction, ghost-cell handling) under `src/operators/`.
   These are oracles retained while their math is ported into catalog
   rules; they are progressively retired as the corresponding rules land
   (see the retirement log in `src/EarthSciDiscretizations.jl`).
   `fv_laplacian` was retired in esd-ecq; its math now lives entirely in
   `discretizations/finite_difference/covariant_laplacian_cubed_sphere.json`.
   New math should land as rules, not operators.

Per the single-pathway rule in [`AGENTS.md`](AGENTS.md), ESD is a
discretization catalog over ESS, not a parallel runtime: no binding
carries a rule evaluator, and every simulation flows through the ESS
`discretize → ArrayOp → eval` pipeline.

See [`docs/REPO_LAYOUT.md`](docs/REPO_LAYOUT.md) for the repository convention
and [`discretizations/README.md`](discretizations/README.md) for how to add a
rule file.

## Development

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
