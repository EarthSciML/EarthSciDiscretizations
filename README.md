# EarthSciDiscretizations

Discretization operators and authoritative discretization rule files for the
EarthSci model stack.

This repository has two roles:

1. A Julia package implementing discretization primitives (finite-volume
   operators on cubed-sphere grids — `fv_laplacian`, `transport_2d`,
   `flux_1d` — PPM transport and reconstruction, ghost-cell handling, and
   grid runtimes for cartesian, lat-lon, vertical, Arakawa-staggered,
   cubed-sphere, MPAS, and DUO topologies under `src/grids/`).
2. A catalog of declarative **discretization rule JSON files** under
   [`discretizations/`](discretizations/), validated against the
   EarthSciSerialization §7 discretization schema and executed by the ESS
   rule engine, with stencil-form lowering (`src/stencil_lowering.jl`) and
   coefficient evaluation (`src/rule_eval.jl`) on the package side.

See [`docs/REPO_LAYOUT.md`](docs/REPO_LAYOUT.md) for the repository convention
and [`discretizations/README.md`](discretizations/README.md) for how to add a
rule file.

## Development

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
