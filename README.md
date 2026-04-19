# EarthSciDiscretizations

Discretization operators and authoritative discretization rule files for the
EarthSci model stack.

This repository has two roles:

1. A Julia package implementing discretization primitives (finite-volume
   operators on cubed-sphere grids, PPM transport, ghost-cell handling,
   the `FVCubedSphere` discretization pipeline, etc.).
2. A catalog of declarative **discretization rule JSON files** under
   [`discretizations/`](discretizations/), validated against the
   EarthSciSerialization §7 discretization schema and executed by the ESS
   rule engine.

See [`docs/REPO_LAYOUT.md`](docs/REPO_LAYOUT.md) for the repository convention
and [`discretizations/README.md`](discretizations/README.md) for how to add a
rule file.

## Development

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```
