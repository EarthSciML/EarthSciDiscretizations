# Covariant-FV (latlon) Path-A production wiring — VERDICT

**Bead:** `esd-6g4.10` (G11: covariant-FV (latlon) production wiring + integration)
**Author:** polecat `gastown.slit` (verdict); polecat `gastown.nux` (resolution)
**Date:** 2026-06-21
**Engine under test:** `EarthSciSerialization` v0.6.0 (`f2b411a6…`, verdict); ESS `origin/main` `78b6a577` (resolution)

## ✅ RESOLUTION — SHIPPED via Option C (ess-gj4)

The mayor adopted **Option C** below. `ess-gj4` ("generic const_array gather
boundary policy") landed on ESS `origin/main` (`78b6a577`) and closed: it adds a
declarative per-const-array, per-dimension boundary policy —
`:periodic` → `mod1(i,N)` wrap, `:clamp` → edge-extend, `:error` → throw
(default) — exposed as a `const_array_boundaries` kwarg on `build_evaluator`. With
it the authored covariant **Laplacian** rule runs **unmodified** through the real
engine; no rule was reformulated and no imperative operator was written.

Path-A production wiring landed in `src/ode_problem.jl` (esd-6g4.10):

1. **Routing.** A latlon GDD that declares `discretizations` (the covariant rules)
   routes through Path A (`_gdd_has_discretizations`); a bare-grid latlon GDD still
   takes Path B. Path selection stays per-family-consistent.
2. **Blocker 1 (rank).** `_grid_primitive_arrays(::LatLonGrid)` now binds every
   metric factor as a `(nlat, nlon)` matrix `M[lat, lon]` (regular grids only),
   gathered `index(name, lat, lon)`. `_unstructured_grid_const_arrays` includes
   `LatLonGrid`; `_inject_grids!` constructs the grid into `loaded`.
3. **Blocker 2 (boundary).** `_grid_const_array_boundaries` emits
   `(:clamp, :periodic)` (lat pole edge-extend; lon periodic wrap) for the metric
   arrays, threaded into `build_evaluator(...; const_array_boundaries=...)`.

**Empirical (real engine, `build_ode_problem` → `build_evaluator`):** build no
longer throws `E_TREEWALK_CONSTARRAY_OOB` at any pole/seam cell; the discrete
Laplace-Beltrami applied to `sin(lat)` converges O(h²) to `−2 sin(lat)/R²` on the
interior — L∞ = 0.0153 / 0.00397 / 0.00100 at N = 16/32/64 (ratios 3.86, 3.96).
Unit coverage: `test/test_ode_problem_covariant_latlon.jl`. Integration solve:
`discretizations/finite_volume/covariant_fv_laplacian_latlon/fixtures/integration/`
(`mms_kind = sin_lat_latlon_laplacian_interior`). The esd-zk9.1 "DECLARATIVE-FEASIBLE"
verdict is upheld **with** Option C engine support (not in the as-authored form
without it — see below). Pole rows carry the documented zero-ghost vs sentinel→self
divergence (§6.2) and are excluded from interior norms.

---

## Verdict (as-authored, pre-ess-gj4): ❌ INFEASIBLE AS AUTHORED — escalated for re-scope

The esd-zk9.2 covariant **Laplacian** rule
(`discretizations/finite_volume/covariant_fv_laplacian_latlon.json`) **cannot be
run through the production ESS engine** (`discretize` → `build_evaluator`). This
is **not** a `build_ode_problem` wiring defect — it is a fundamental
incompatibility between the rule's formulation and the ESS engine's `const_array`
gather semantics, proven empirically below. Per the bead's **DECLARATIVE-OR-FAIL**
clause ("If it CANNOT be expressed declaratively … STOP and report the precise gap
for escalation"), no production wiring was landed and no imperative/hand-coded
operator was written.

The covariant **gradient** rule (`covariant_fv_gradient_latlon.json`) **does** run
through the engine once the metric binding is corrected (see Blocker 1) — it
gathers the metric at the cell centre only.

## TL;DR

| Rule | Through real ESS engine? | Why |
|---|---|---|
| `covariant_fv_gradient_latlon_{t1,t2}` | ✅ works (after 2-D binding fix) | metric gathered at `[lat,lon]` (centre) only |
| `covariant_fv_laplacian_latlon` | ❌ `E_TREEWALK_CONSTARRAY_OOB` | connection terms gather metric at `lat±1`/`lon±1` offsets; ESS never wraps/ghosts `const_array` gathers |

Two independent blockers, both in the host↔engine metric binding (not the
rule's *expressibility* per se, but its *runnability* over the existing engine):

1. **Rank mismatch (`E_TREEWALK_CONSTARRAY_NDIM`).** `_grid_primitive_arrays(::LatLonGrid)`
   (`src/ode_problem.jl:862-883`) binds every metric array as a **1-D vector** of
   length `nc = nlon*nlat`, but the rules gather them with **two** indices
   `index(g_xx, lat, lon)`. ESS requires `length(index args) == ndims(array)` and
   throws otherwise. (Also: `_unstructured_grid_const_arrays`,
   `src/ode_problem.jl:805-812`, only emits arrays for `DuoGrid`/`MpasGrid` —
   `LatLonGrid` is excluded at line 808, so today the binding is never even reached.)
   *Fixable* host-side: bind as `(nlat, nlon)` matrices. Necessary for both rules.

2. **Boundary out-of-bounds (`E_TREEWALK_CONSTARRAY_OOB`) — the fatal one.** The
   Laplacian's four connection-correction terms (rule lines 56-135) form centered
   differences of the metric: `(Jg_xx[lon+1] − Jg_xx[lon-1])/(2·dlon)`,
   `(Jg_yy[lat+1] − Jg_yy[lat-1])/(2·dlat)`, and the two `Jg_xe` cross terms. These
   gather the metric `const_arrays` at **offset** indices. ESS applies **no periodic
   wrap and no ghost-cell substitution to `const_array` gathers** — an out-of-range
   `const_array` index is a hard throw. So at every boundary cell (the periodic-lon
   edge `lon-1` at `lon=1`, and the non-periodic-lat poles `lat±1`) the build aborts.
   This is **not** host-fixable without changing the rule or the engine.

## Background: the two discretization paths

`build_ode_problem` (`src/ode_problem.jl:57`) chooses between:

- **Path A** — the declarative rule engine: `EarthSciSerialization.discretize(esm; lift_1d_arrayop=true)`
  then `build_evaluator(disc; const_arrays=…)` (`src/ode_problem.jl:84,100`). This
  is where the esd-zk9 declarative rules would run.
- **Path B** — PDESystem/MTK: `_build_path_b` (`src/ode_problem.jl:122`) flattens
  to a `ModelingToolkit.PDESystem` and calls the curvilinear FD chain-rule pipeline
  `discretize(sys, grid; xi_axis, eta_axis)`. **This path ignores the Path-A
  declarative rules entirely.**

Today `latlon` is forced onto Path B by `_CURVILINEAR_FAMILIES = {"latlon"}`
(`src/ode_problem.jl:7`) via `_has_curvilinear_domain` →
`return _build_path_b(...)` (`src/ode_problem.jl:71-72`). G11's goal was to remove
that and route latlon covariant through Path A so the declarative rules execute in
production. The rule injection already translates the authored geographic axis
names to canonical loop indices — `_inject_rules!` maps `lat→i`, `lon→j` for
`grid_family == "latlon"` (`src/ode_problem.jl:655-660`).

## Root cause: ESS `const_array` gather semantics

The decisive asymmetry is in the ESS tree-walk gather resolver
`_resolve_indices` (`EarthSciSerialization/src/tree_walk.jl`):

- **State-variable gathers** `index(u, …)` (`tree_walk.jl:2163-2181`): rank comes
  from discovered cell bounds `(lo,hi)`; an out-of-range index returns
  `NumExpr(0.0)` — the **ghost-cell convention** (`tree_walk.jl:2174-2175`).
  Periodic dimensions are additionally folded by `_apply_periodic_folding!`
  (`discretize.jl:597-648`) — but **only for variables carrying `shape`+`grid`
  metadata** (`discretize.jl:546-563`).
- **`const_array` gathers** `index(Jg_xx, …)` (`tree_walk.jl:2185-2200`): rank
  comes from `ndims(array)`; `length(index args) == ndims` is enforced
  (`E_TREEWALK_CONSTARRAY_NDIM`), and any out-of-range index throws
  `E_TREEWALK_CONSTARRAY_OOB` (`tree_walk.jl:2194-2198`). **No ghost, no clamp,
  no periodic wrap — ever** (the periodic folder excludes `const_arrays` because
  they are not `shape`+`grid` variables).

So a rule that reads the metric at offset indices across a grid boundary is
unrunnable: state offsets are ghosted/periodic-folded, but the metric `const_array`
offsets are not. The gradient rule reads the metric only at the centre, so it is
unaffected.

This contradicts the load-bearing assumption in the esd-zk9.1
`ORACLE_CHARACTERIZATION.md` ("ESS gathers them via `index(name, lat, lon)` and
const-folds … flat per-cell layout lon-fastest"), and the matching comment at
`src/ode_problem.jl:858-861`. ESS does **not** flatten a multi-index into a 1-D
`const_array` (it errors), and it does **not** wrap/ghost `const_array` offset
reads at boundaries. The esd-zk9 rules were validated only by the **test-local
interpreter** in `test/test_covariant_fv_metric_conformance.jl`, which resolves
gathers itself (manual flatten + sentinel→self boundary) and therefore never
exercised these two engine behaviors. The bead anticipated exactly this risk
("currently only UNIT via a test-local interpreter").

## Empirical evidence

Both spikes drive the **authored** rule JSON through the **real** ESS engine,
replicating the production injection (`_subst_axis_names` lat→i/lon→j), and bind
the metric as corrected 2-D `(nlat,nlon)` matrices so Blocker 1 is bypassed and
Blocker 2 is observable. `nlat=nlon=6`, `R=1`, lon periodic, lat non-periodic.

### Laplacian — FAILS

```
=== STEP 1: discretize (rule engine) ===
discretize OK
=== STEP 2: build_evaluator (per-cell expansion + const gathers) ===
build_evaluator THREW: EarthSciSerialization.TreeWalkError
E_TREEWALK_CONSTARRAY_OOB: const array 'Jg_xe' index 0 out of range 1..6 in dim 2
```

`Jg_xe`, dim 2 (lon), index 0 = `index(Jg_xe, lat, lon-1)` at `lon=1`. lon is
periodic, but the `const_array` gather is not wrapped → hard throw.

### Gradient (t1) — WORKS

```
=== GRADIENT t1: discretize ===
discretize OK
=== GRADIENT t1: build_evaluator ===
build OK; len(u0)=36 finite=true
GRADIENT SPIKE: RUNS THROUGH REAL ESS ENGINE ✅
```

Reproduction scripts are reproduced verbatim in the Appendix.

## Why this is not a `build_ode_problem` wiring fix

The only host-side levers are: which path latlon takes (`_CURVILINEAR_FAMILIES`),
constructing a `LatLonGrid` in `_inject_grids!`, and the shape/contents of the
metric `const_arrays`. None of these can give a `const_array` gather a periodic
or ghost boundary — that behavior lives entirely inside the ESS engine and is
reserved for `shape`+`grid` state variables. Correcting the binding to 2-D
(Blocker 1) is necessary but not sufficient; the Laplacian still aborts on
Blocker 2 at every edge cell. Routing only the gradient to Path A is not viable
either, because path selection is per-GDD-**family**, not per-operator — a latlon
GDD that uses `laplacian` would crash.

## Resolution options (for mayor decision)

**A. Reformulate the Laplacian as precomputed centre-gathered 9-point weights**
(host data; declarative; existing vocabulary). Emit nine per-cell weight
`const_arrays` `W_C,W_N,W_S,W_E,W_W,W_NE,W_NW,W_SE,W_SW` (the exact stencil weights
the oracle `grid_assembly.jl` already produces — esd-zk9.1 §5 verified all nine
columns to 0.0 error), and author the rule as
`Σ_k index(W_k, lat, lon) · index($u, lat+δ_k, lon+δ_k)`. All metric/weight gathers
are at the centre (no OOB); the only offsets are on `$u` (ghosted/periodic-folded
by the engine). Boundary policy is baked into the weights by the host. *This
re-opens the esd-zk9.2 rule deliverable and requires re-verifying the new rule
against the oracle goldens; it is not a wiring task.*

**B. Reformulate the connection terms with centre-gathered metric-derivative
arrays.** Replace `(Jg_xx[lon+1]−Jg_xx[lon-1])/(2·dlon)` etc. with a single
centre gather `index(dJg_xx_dlon, lat, lon)` of a host-precomputed derivative
array (the grid already exposes `metric_dgij_dxk`, `src/grids/latlon.jl:768`).
Same scope caveat as A.

**C. Add GENERIC ESS engine support for `const_array` boundary handling** —
periodic wrap and/or ghost for `const_array` gathers, mirroring the state-variable
path. This is the only option that runs the **authored** rules unmodified, but it
is an `EarthSciSerialization` change (separate package), and ghost-0 for a metric
at a non-periodic pole is physically wrong (periodic-lon wrap would be correct).
The bead permits "GENERIC engine support," but this is cross-package and needs
ESS-side design.

## Recommendation

Pursue **Option A** (precomputed centre-gathered weights) as a re-scoped
rule-authoring bead on the covariant-FV thread (sibling of esd-zk9.2), since the
oracle already produces the exact weights and the form is provably engine-runnable
(it reduces the Laplacian to the same centre-coefficient × state-stencil shape the
gradient already uses successfully). The corrected **2-D metric binding** (Blocker
1) and the **gradient Path-A wiring** are small and can ride along once the
Laplacian form is settled, so that path selection stays per-family-consistent.

Either way, the esd-zk9.1 `ORACLE_CHARACTERIZATION.md` "DECLARATIVE-FEASIBLE"
verdict should be amended: it is feasible, but **not in the metric-arrays + in-rule
finite-difference form esd-zk9.2 authored** — only in a centre-gathered form
(A/B), or with new ESS engine support (C).

## Appendix — reproduction

Env: `julia --project=.` in this worktree; ESS v0.6.0 (`f2b411a6…`).

### Spike 1 — Laplacian (fails)

```julia
import EarthSciDiscretizations as ESD
import EarthSciSerialization as ESS
import JSON
REPO = @__DIR__  # point at repo root in practice
nlat, nlon, R = 6, 6, 1.0
dlat, dlon = pi/nlat, 2pi/nlon
spec = JSON.parsefile(joinpath(REPO,"discretizations/finite_volume/covariant_fv_laplacian_latlon.json"))["discretizations"]["covariant_fv_laplacian_latlon"]
replacement = ESD._subst_axis_names(spec["replacement"], Dict("lat"=>"i","lon"=>"j"))
esm = Dict{String,Any}(
  "esm"=>"0.4.0","metadata"=>Dict("name"=>"spike"),
  "grids"=>Dict("domain"=>Dict("family"=>"latlon","dimensions"=>Any[
     Dict("name"=>"lat","size"=>nlat,"periodic"=>false,"spacing"=>"uniform"),
     Dict("name"=>"lon","size"=>nlon,"periodic"=>true,"spacing"=>"uniform")])),
  "models"=>Dict("M"=>Dict("grid"=>"domain","variables"=>Dict(
     "u"=>Dict("type"=>"state","shape"=>["lat","lon"],"location"=>"cell_center","default"=>0.0),
     "dlat"=>Dict("type"=>"parameter","default"=>dlat),
     "dlon"=>Dict("type"=>"parameter","default"=>dlon)),
   "equations"=>Any[Dict("lhs"=>Dict("op"=>"D","args"=>["u"],"wrt"=>"t"),
                         "rhs"=>Dict("op"=>"laplacian","args"=>["u"]))])),
  "rules"=>Any[Dict("name"=>"covariant_fv_laplacian_latlon","pattern"=>spec["applies_to"],"replacement"=>replacement)])
grid = ESD._latlon(nlon=nlon,nlat=nlat,R=R)
ginv, jac = ESD.metric_ginv(grid), ESD.metric_jacobian(grid)
nc = ESD.n_cells(grid)
flat2mat(v) = permutedims(reshape(Float64.(v), nlon, nlat),(2,1))  # (nlat,nlon)
ca = Dict{String,AbstractArray{Float64}}(
  "g_xx"=>flat2mat([ginv[k,1,1] for k in 1:nc]), "g_yy"=>flat2mat([ginv[k,2,2] for k in 1:nc]),
  "g_xe"=>flat2mat([ginv[k,1,2] for k in 1:nc]), "invJ"=>flat2mat([1/jac[k] for k in 1:nc]),
  "Jg_xx"=>flat2mat([jac[k]*ginv[k,1,1] for k in 1:nc]), "Jg_yy"=>flat2mat([jac[k]*ginv[k,2,2] for k in 1:nc]),
  "Jg_xe"=>flat2mat([jac[k]*ginv[k,1,2] for k in 1:nc]))
disc = ESS.discretize(esm; lift_1d_arrayop=true)          # OK
ESS.build_evaluator(disc; initial_conditions=Dict{String,Float64}(), const_arrays=ca)  # THROWS E_TREEWALK_CONSTARRAY_OOB
```

### Spike 2 — Gradient t1 (works)

Identical setup with the gradient rule
(`covariant_fv_gradient_latlon_t1`, `rhs = grad(u, dim=t1)`) and the two centre
metric arrays `dxi_dt1`, `deta_dt1` from `coord_jacobian(grid,:lon_lat)`:
`build_evaluator` returns a finite RHS over all 36 cells.
