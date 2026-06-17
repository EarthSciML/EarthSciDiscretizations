# EarthSciDiscretizations.jl

**EarthSciDiscretizations.jl (ESD)** is a *declarative catalog of PDE-discretization rules* over [EarthSciSerialization (ESS)](https://github.com/EarthSciML/EarthSciSerialization.jl), plus a multi-family grid runtime and a top-level `build_ode_problem` entry point that turns a PDE component (`.esm`) into a ready-to-`solve` `ODEProblem`.

ESD is three things, in priority order:

1. **A catalog of declarative discretization-rule JSON files** under `discretizations/` (finite-difference and finite-volume stencils, reconstructions, flux forms, and boundary conditions). The rules are applied by the **ESS rule engine** via the canonical `discretize → ArrayOp → eval` pipeline; ESD carries no rule evaluator of its own.
2. **A multi-family grid runtime** (`src/grids/`): `CartesianGrid`, `VerticalGrid`, `LatLonGrid`, `ArakawaGrid` (A/B/C/D/E staggers), `CubedSphereGrid`, `MpasGrid`, and `DuoGrid`, mirrored by Python/Rust/TypeScript bindings (see [`GRIDS_API.md`](https://github.com/EarthSciML/EarthSciDiscretizations.jl/blob/main/docs/GRIDS_API.md)).
3. **A top-level ESM → ODEProblem constructor**, `build_ode_problem` (see below), which loads a PDE `.esm`, optionally merges a Grid Discretization Descriptor (`*.gdd.json`), runs the ESS pipeline, and returns `(prob::ODEProblem, var_map::Dict{String,Int})`.

The mental model is: **declarative PDE (`.esm`) + grid descriptor (`*.gdd.json`) → `build_ode_problem` → `ODEProblem` → `solve`.**

The finite-volume cubed-sphere transport algorithms in the catalog are based on the Lin-Rood (1996) dimensionally-split scheme with Piecewise Parabolic Method (PPM) reconstruction (Colella & Woodward 1984), following Putman & Lin (2007); vertical remapping follows Lin (2004). See [Finite-Volume Method](@ref) for full references.

## Features

- **Declarative discretization rules** authored as closed `arrayop`/`makearray` ASTs in the §4.2 ESS op vocabulary (no scheme-specific Julia kernels). The catalog covers centered finite differences, upwind, flux limiters, PPM reconstruction, Lin-Rood transport, and grid-family-specific Laplacians. See [Operators](@ref) and [Finite-Volume Method](@ref).
- **Seven grid families** with a keyword-constructor submodule `EarthSciDiscretizations.grids.<family>(; …)` (plus `CubedSphereGrid(Nc; …)`), automatic metric-tensor computation, ghost-cell management, and Arakawa staggering.
- **`build_ode_problem` entry point** that routes a PDE through the ESS rule engine (Path A: cartesian/vertical/arakawa/mpas/duo) or the ESS `PDESystem` pipeline (Path B: latlon/cubed_sphere) and returns a `SciMLBase.ODEProblem`.
- **Declarative boundary and initial conditions** carried on the model and applied by the ESS rule engine; expression initial conditions are rewritten to ESS IC-arrayop equations.
- **Cross-binding grid API** (Julia/Python/Rust/TypeScript) defined normatively in [`GRIDS_API.md`](https://github.com/EarthSciML/EarthSciDiscretizations.jl/blob/main/docs/GRIDS_API.md).

## Quick Start

`build_ode_problem` is the primary entry point: give it a PDE `.esm` and a grid
descriptor and it returns an `ODEProblem` you can `solve`.

```@example quickstart
using EarthSciDiscretizations
using OrdinaryDiffEqDefault: solve

repo = dirname(dirname(pathof(EarthSciDiscretizations)))
esm  = joinpath(repo, "test", "fixtures", "diffusion_1d.esm")
gdd  = joinpath(repo, "discretizations", "gdd", "cartesian_1d_periodic_n16.gdd.json")

# .esm + GDD → ODEProblem + state-name → index map
prob, var_map = build_ode_problem(esm; grid_ref = gdd)

# Seed a nontrivial initial condition (16 cells, periodic [0,1]).
N  = 16
dx = 1.0 / N
for i in 1:N
    prob.u0[var_map["u[$i]"]] = sinpi(2 * (i - 0.5) * dx)
end

sol = solve(prob; saveat = [0.0, 0.05])
println("retcode:  ", sol.retcode)
println("u[1] @ t=0:    ", round(sol.u[1][var_map["u[1]"]];   digits = 6))
println("u[1] @ t=0.05: ", round(sol.u[end][var_map["u[1]"]]; digits = 6))
```

For a full walkthrough — including how to write the `.esm`, what a GDD
contains, Path A vs Path B routing, and curvilinear (`remake`) cases — see
[Getting started: solve a PDE](@ref).

## The `build_ode_problem` entry point

```julia
build_ode_problem(esm_path; grid_ref="", reader_fn=nothing, extra_ics=Dict()) -> (prob, var_map)
```

Loads the PDE component `.esm`, optionally merges the GDD at `grid_ref`, runs
the ESS discretization pipeline, and returns a `SciMLBase.ODEProblem` together
with `var_map::Dict{String,Int}` mapping each scalar state name (e.g.
`"u[1]"`) to its index in `prob.u0`. The constructor never invokes a solver —
you `solve` it yourself, with `OrdinaryDiffEqDefault` as the intended solver
dependency. See [Getting started: solve a PDE](@ref) for the kwargs and the
two routing paths.

## Package Contents

- [Getting started: solve a PDE](@ref) -- end-to-end `.esm` + GDD → `build_ode_problem` → `solve` walkthrough (start here)
- [Finite-Volume Method](@ref) -- mathematical foundations and how a rule encodes an FV operator
- [Operators](@ref) -- the closed §4.2 op vocabulary that rule replacements use
- [Tutorial: Authoring a rule](@ref) -- end-to-end walkthrough for writing a new discretization rule in the closed-AST pattern

!!! note "Two documentation sites"
    This Documenter site covers the **Julia package API** (`build_ode_problem`,
    grids, the operator vocabulary, and how to author a rule). A companion
    **[catalog browser](https://EarthSciML.github.io/EarthSciDiscretizations.jl/catalog/)**
    (auto-generated from `discretizations/` and `src/grids/`) renders a page
    per rule and per grid family with stencil and convergence plots. Use this
    site for "how do I use ESD"; use the catalog browser for "what rules and
    grids exist."
