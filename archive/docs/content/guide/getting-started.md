---
title: "Getting started: solve a PDE"
slug: "getting-started"
description: "From a PDE written as an .esm file to a solved ODEProblem: the build_ode_problem workflow exercised throughout the ESD test suite."
---

This page takes you from a partial differential equation, written as an
`.esm` file, to a solved `ODEProblem`. It is the workflow exercised
throughout the test suite (`test/test_build_ode_problem.jl`,
`test/integration_cases/*`).

The whole pipeline is one call plus a `solve`:

> **declarative PDE (`.esm`) + grid descriptor (`*.gdd.json`)
> → `build_ode_problem` → `ODEProblem` → `solve`**

ESD constructs the spatially-discretized `ODEProblem`; you choose the
time integrator. ESD itself depends on no solver — add one to your own
project; `OrdinaryDiffEqDefault` is the recommended default.

## A complete example (Path A)

```julia
using EarthSciDiscretizations
using OrdinaryDiffEqDefault: solve

repo = dirname(dirname(pathof(EarthSciDiscretizations)))
esm  = joinpath(repo, "test", "fixtures", "diffusion_1d.esm")
gdd  = joinpath(repo, "discretizations", "gdd", "cartesian_1d_periodic_n16.gdd.json")

prob, var_map = build_ode_problem(esm; grid_ref = gdd)

# var_map maps each scalar state name to its index in prob.u0.
N  = 16
dx = 1.0 / N
for i in 1:N
    prob.u0[var_map["u[$i]"]] = sinpi(2 * (i - 0.5) * dx)
end

sol = solve(prob; saveat = [0.0, 0.05])
println("retcode: ", sol.retcode)
println("u[1] @ t=0:    ", round(sol.u[1][var_map["u[1]"]];   digits = 6))
println("u[1] @ t=0.05: ", round(sol.u[end][var_map["u[1]"]]; digits = 6))
```

That solves the 1D diffusion equation `∂u/∂t = ∂²u/∂x²` on a periodic
domain of 16 cells. Each step below explains one piece of it.

## The pieces

### 1. The PDE — an `.esm` file

`diffusion_1d.esm` is a complete PDE component: a model with one state
variable `u`, a spacing parameter `dx`, and one equation
`∂u/∂t = d2(u, dim=x)` (the continuous operators are written as an AST):

```json
{
  "esm": "0.4.0",
  "metadata": {"name": "diffusion_1d_pde"},
  "models": {
    "M": {
      "grid": "domain",
      "variables": {
        "u":  {"type": "state",     "shape": ["x"], "location": "cell_center", "default": 0.0},
        "dx": {"type": "parameter"}
      },
      "equations": [
        { "lhs": {"op": "D",  "args": ["u"], "wrt": "t"},
          "rhs": {"op": "d2", "args": ["u"], "dim": "x"} }
      ]
    }
  }
}
```

The `.esm` is *grid-agnostic*: it names a grid (`"domain"`) but does not
say how big it is or which discretization rules apply. That is the GDD's
job.

### 2. The grid + rules — a Grid Discretization Descriptor (`*.gdd.json`)

The GDD is the user's grid-and-rule selector. It declares the concrete
grid (extent, spacing, boundary conditions) and the catalog rules to
apply. `cartesian_1d_periodic_n16.gdd.json`:

```json
{
  "esm": "0.5.0",
  "kind": "grid_discretization_descriptor",
  "grids": {
    "domain": {
      "spatial": { "x": { "min": 0.0, "max": 1.0, "grid_spacing": 0.0625 } },
      "boundary_conditions": [ { "type": "periodic", "dimensions": ["x"] } ]
    }
  },
  "discretizations": {
    "centered_2nd_uniform":       { "ref": "../finite_difference/centered_2nd_uniform.json" },
    "centered_2nd_deriv_uniform": { "ref": "../finite_difference/centered_2nd_deriv_uniform.json" },
    "upwind_1st":                 { "ref": "../finite_difference/upwind_1st.json" }
  }
}
```

`build_ode_problem` merges this into the `.esm` (`grid_ref` may be an
absolute path or one relative to the `.esm`'s directory): the `grids`
block sizes the named grid and injects per-cell parameters (e.g. `dx`),
and each `discretizations` entry is a `ref` to a catalog rule JSON under
`discretizations/`. The rules are *not* a registry lookup — the GDD
points at the rule files by relative path.

### 3. The call — `build_ode_problem`

```julia
build_ode_problem(esm_path; grid_ref="", reader_fn=nothing, extra_ics=Dict()) -> (prob, var_map)
```

| Keyword | Purpose |
|---|---|
| `grid_ref` | Path to the GDD to merge (absolute, or relative to `esm_path`). Empty string ⇒ no grid merge (the `.esm` must already be sized). |
| `reader_fn` | Loader callback for file-backed grids (MPAS `.nc`/mesh). Used by the MPAS path; `nothing` otherwise. |
| `extra_ics` | `Dict{String,Float64}` of extra numeric per-cell initial values merged into the IC set (e.g. field-pipeline column data). |

It runs the ESS discretization pipeline and returns:

- `prob::SciMLBase.ODEProblem` — the spatially-discretized system. **No
  solver is invoked inside the constructor.**
- `var_map::Dict{String,Int}` — maps each scalar state name to its index
  in `prob.u0`.

### 4. Indexing the state — `var_map`

`var_map` is how you read and write individual cells of `prob.u0`. For a
state `u` with `N` cells the keys are `"u[1]" … "u[N]"`:

```julia
prob.u0[var_map["u[1]"]]       # value at cell 1
prob.u0[var_map["u[$i]"]] = …  # set cell i
```

After solving, `sol.u[k]` is the full state vector at save point `k`, so
`sol.u[end][var_map["u[$i]"]]` is the final value at cell `i`. (For a
multi-variable model the keys are `"<name>[<cell>]"` per variable.)

### 5. The solve

ESD hands back a standard SciML `ODEProblem`, so any
`OrdinaryDiffEq`/SciML solver works. `OrdinaryDiffEqDefault.solve`
auto-selects an integrator and is the recommended default (add it to
your own project — ESD does not depend on a solver):

```julia
using OrdinaryDiffEqDefault: solve
sol = solve(prob; saveat = [0.0, 0.05])
```

## Path A vs Path B (how the grid family routes)

`build_ode_problem` routes on the grid family declared in the GDD:

| Path | Grid families | Pipeline | IC handling |
|---|---|---|---|
| **A** (default) | `cartesian`, `vertical`, `arakawa`, `mpas`, `duo`, and the no-grid case | ESS rule engine: `discretize → ArrayOp → build_evaluator` | ICs are carried into `prob.u0` (see below) |
| **B** (curvilinear) | `latlon` | ESS `PDESystem` pipeline: `load → flatten → PDESystem → discretize(sys, grid)` | Produces a **zero IC** — inject your physical IC with `remake` (see below) |

The example above is Path A. For a curvilinear family you set
`family: "latlon"` in the GDD; the example GDD is
`test/fixtures/curvilinear/latlon_diffusion.gdd.json`, e.g.:

```json
{ "grids": { "sphere": { "family": "latlon", "R": 1.0,
    "spatial": {
      "lon": {"min": -3.14159, "max": 3.14159, "grid_spacing": 1.0472},
      "lat": {"min": -1.5708,  "max": 1.5708,  "grid_spacing": 0.5236}
    } } } }
```

`build_ode_problem` on that GDD returns a solvable `ODEProblem` whose
state vector has one entry per cell (here `nlon=6 × nlat=6 = 36`), as
exercised in `test/test_ode_problem_curvilinear.jl`.

### Path B and the initial condition

The Path-B `PDESystem` pipeline produces a problem whose initial
condition is **zero** (the spatial BCs are baked into the discretized
operator, not the IC vector). To run it with a physical initial field,
build your IC vector over the grid's cell centres and swap it in with
`SciMLBase.remake`, which preserves MTK's initialization data while
replacing `u0`:

```julia
import SciMLBase
prob_template, var_map = build_ode_problem(esm; grid_ref = gdd)  # zero IC
u0 = my_initial_condition(...)                                   # length == length(prob_template.u0)
prob = SciMLBase.remake(prob_template; u0 = u0)                  # inject the IC
sol  = solve(prob; reltol = 1e-6, abstol = 1e-8, save_everystep = false)
```

## Initial conditions

There are two kinds of initial condition on a Path-A model:

- **Expression ICs** — a model whose `initial_conditions.type ==
  "expression"` gives each field a formula over the spatial coordinates,
  e.g. `u(x,0) = sin(2π x)` written as an AST. `build_ode_problem`
  *rewrites* these into the model's `initialization_equations` and lets
  the **ESS IC-arrayop engine** evaluate them, binding `coord_<dim>`
  arrays of cell centres. (They are not evaluated by ESD itself.)

  ```json
  "initial_conditions": {
    "type": "expression",
    "values": { "u": {"op": "sin", "args": [{"op": "*", "args": [6.283185307179586, "x"]}]} }
  }
  ```

- **Numeric per-cell ICs** — concrete per-cell values (e.g. `dz[k]`
  column widths injected by the GDD, or anything you pass via
  `extra_ics`) are placed directly into `prob.u0`.

## Boundary conditions and current limits

Boundary conditions are declared on the model's `boundary_conditions`
block (or, for whole-axis conditions like periodicity, in the GDD's
`grids.*.boundary_conditions`). The ESS rule engine applies them. A
per-side block looks like:

```json
"boundary_conditions": {
  "u_xmin": {"variable": "u", "kind": "neumann", "side": "xmin", "value": 1.0},
  "u_xmax": {"variable": "u", "kind": "neumann", "side": "xmax", "value": 0.0}
}
```

<div class="callout callout-note">
<strong>Which BC kinds work today.</strong>
As of the current trunk, <strong>periodic</strong>, <strong>zero-value
Dirichlet</strong>, <strong>nonzero-Neumann</strong>, and
<strong>Robin</strong> BCs all solve end-to-end through
<code>build_ode_problem</code> (see the <code>bc_ic</code> goldens in
<code>test/test_bc_ic_goldens.jl</code>). The nonzero-Neumann and Robin
ghost rules
(<code>discretizations/finite_difference/{neumann_bc,robin_bc}.json</code>)
are wired into the integration pipeline by <code>_inject_bc_rules!</code>;
the ESS rule engine fires them on the synthetic <code>bc</code> node and the
makearray-region BC lowering (ess-hjg) splices the rewritten ghost into the
boundary regions (esd-6k1).
</div>

## Integral / PIDE terms

ESD can discretize partial integro-differential equations. An `integral`
operator in the `.esm` RHS is lowered by the pipeline to a quadrature
over the named axis:

```json
{ "op": "integral", "args": ["u"], "var": "x",
  "lower": {"op": "const", "value": 0.0},
  "upper": {"op": "const", "value": 1.0} }
```

For example `∂u/∂t = -∫₀¹ u dx` (fixture `test/fixtures/pide_integral_1d.esm`)
builds and solves like any other Path-A model; the midpoint rule is exact
for linear integrands. See the PIDE test in
`test/test_build_ode_problem.jl`.

## Where to go next

- [The finite-volume method]({{< ref "/guide/finite-volume-method" >}}) —
  how a rule's pattern match + closed `arrayop` replacement encode an
  operator.
- [Operators]({{< ref "/guide/operators" >}}) — the closed §4.2 op
  vocabulary that every rule replacement uses.
- [Authoring a rule]({{< ref "/guide/authoring-a-rule" >}}) — write your
  own discretization rule.
- [`GRIDS_API.md`]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md) — the
  grid families and their constructor options.
