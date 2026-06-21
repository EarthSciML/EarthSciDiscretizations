---
title: "Guide"
description: "How to use EarthSciDiscretizations — from a PDE to a solved ODEProblem, the operator vocabulary, the finite-volume method, and authoring a rule."
---

This is the **user and authoring guide** for
EarthSciDiscretizations.jl (ESD). The companion
[catalog browser]({{< ref "/" >}}) answers *what rules and grids exist*;
this guide answers *how to use ESD* and *how a rule works*.

ESD is three things, in priority order:

1. **A catalog of declarative discretization-rule JSON files** under
   [`discretizations/`]({{< param repoURL >}}/tree/main/discretizations)
   (finite-difference and finite-volume stencils, reconstructions, flux
   forms, and boundary conditions). The rules are applied by the
   **EarthSciSerialization (ESS) rule engine** via the canonical
   `discretize → ArrayOp → eval` pipeline; ESD carries no rule evaluator
   of its own.
2. **A multi-family grid runtime** (`src/grids/`): `CartesianGrid`,
   `VerticalGrid`, `LatLonGrid`, `ArakawaGrid` (A/B/C/D/E staggers),
   `MpasGrid`, and `DuoGrid`, mirrored by Python/Rust/TypeScript bindings
   (see [`GRIDS_API.md`]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md)).
3. **A top-level ESM → ODEProblem constructor**, `build_ode_problem`,
   which loads a PDE `.esm`, optionally merges a Grid Discretization
   Descriptor (`*.gdd.json`), runs the ESS pipeline, and returns
   `(prob::ODEProblem, var_map::Dict{String,Int})`.

The mental model is: **declarative PDE (`.esm`) + grid descriptor
(`*.gdd.json`) → `build_ode_problem` → `ODEProblem` → `solve`.**

## The `build_ode_problem` entry point

```julia
build_ode_problem(esm_path; grid_ref="", reader_fn=nothing, extra_ics=Dict()) -> (prob, var_map)
```

Loads the PDE component `.esm`, optionally merges the GDD at `grid_ref`,
runs the ESS discretization pipeline, and returns a
`SciMLBase.ODEProblem` together with `var_map::Dict{String,Int}` mapping
each scalar state name (e.g. `"u[1]"`) to its index in `prob.u0`. The
constructor never invokes a solver — you `solve` it yourself, with the
time integrator of your choice. ESD itself takes no solver dependency;
`OrdinaryDiffEqDefault` is the recommended default to add to your own
project.

## Pages

- **[Getting started: solve a PDE]({{< ref "/guide/getting-started" >}})** —
  the end-to-end `.esm` + GDD → `build_ode_problem` → `solve`
  walkthrough. **Start here.**
- **[The finite-volume method]({{< ref "/guide/finite-volume-method" >}})** —
  the mathematical foundations and how a rule's pattern match + closed
  `arrayop` replacement encode an FV operator.
- **[Operators]({{< ref "/guide/operators" >}})** — the closed §4.2 op
  vocabulary that every rule replacement uses.
- **[Authoring a rule]({{< ref "/guide/authoring-a-rule" >}})** — the
  conceptual closed-AST walkthrough for writing a new discretization
  rule. For the full repository-infrastructure version (paths, CI layers,
  walker registration) see the
  [Add a new discretization rule]({{< ref "/tutorials/add-a-rule" >}})
  contributor tutorial.
