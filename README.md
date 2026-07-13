# EarthSciDiscretizations

The standard library of discretization rules for the
[EarthSciAST](https://github.com/EarthSciML/EarthSciAST) (ESS/ESM)
format: grids, finite-difference/finite-volume rewrite rules, regridding and
reprojection expressions — together with the cross-language conformance goldens,
convergence suites, and MMS tests that pin their behavior across every ESS binding.

This is a **pure data repository**. Every library entry is a `.esm` document; there is
no evaluator, rule engine, or numeric kernel here (see [AGENTS.md](AGENTS.md),
"single pathway"). The ESM spec (esm-spec.md §9.6.8) ships no discretization rules
itself and delegates the standard library to this repo.

## The layering

Library files compose by reference (esm-spec §9.7): each layer is a valid ESM
template-library file importing the one below it, and every file is
resolution-generic through load-time **metaparameters**.

```
grids/<grid>/grid.esm             index sets + geometry (metaparameters: N, NLON, ...)
  └─ stencils/<name>.esm          interior operator, match-less template
       └─ rules/<name>.esm        stencil + boundary conditions = the complete
                                  auto-applied rewrite rule on spatial D
            └─ problems/<name>.esm   runnable driver with inline MMS tests
                                     (also: your model, importing the rule)
```

A consuming model needs one import and one `D`:

```json
"expression_template_imports": [
  { "ref": ".../grids/cartesian_uniform_1d/rules/central_D2_zero_grad_bc.esm",
    "bindings": { "N": 128 } }
],
"equations": [
  { "lhs": { "op": "D", "args": ["u"], "wrt": "t" },
    "rhs": { "op": "*", "args": ["alpha",
             { "op": "D", "args": [{ "op": "D", "args": ["u"], "wrt": "x" }], "wrt": "x" }] } }
]
```

Boundary conditions live **inside** the rule (interior region = the stencil aggregate,
boundary-face regions = the BC); there is no separate boundary-condition declaration
anywhere (esm-spec §9.6.8). Rebinding the metaparameters at the import edge is the whole
convergence story — the same files serve every resolution.

The *extent* is yours too. A grid's real-valued geometry is not baked in: the origins
and spacings are consumer-supplied free names (`x0`/`dx` on `cartesian_uniform_1d`,
`lon0_deg`/`dlon_deg`/`lat0_deg`/`dlat_deg` + `R_sphere` on `latlon`) that you define as
ordinary model variables — the MPAS keyed-factor pattern. For a domain `[a, b]` you set
`x0 = a` and `dx = (b − a)/N` (an observed dividing by the metaparameter name `N`, so a
convergence sweep stays consistent). The same rule file that runs on the unit interval
runs on `x ∈ [−1.5, 2.5]` at observed order 2.00
(`problems/heat_1d_zero_grad_nonunit.esm`). And because every rule carries a §9.6.1
`where` shape constraint that travels with §9.7.7 import renaming, you can import one
grid+rule family **twice** in one model — each instance scoped to its own mesh with its
own spacing (`problems/two_cartesian_grids_coexist.esm`).

## Layout

| Directory | Contents |
|---|---|
| `grids/<grid>/` | One directory per grid: `grid.esm` (the constructor), `stencils/`, `rules/` — the same scheme takes a different form on different grids |
| `regridding/` | Conservative/interpolating remap expressions (coupling transforms) |
| `reprojection/` | Coordinate-transform template fragments (lon-lat, Lambert conformal, …) |
| `problems/` | Runnable driver models with inline §6.6.5 MMS tests |
| `tests/conformance/` | `ast/` (byte-identical post-lowering ASTs, all five bindings), `simulation/` (MMS), `convergence/` (observed order), `regridding/`, `reprojection/` |
| `tests/invalid/` | `lint/` (policy-lint fixtures) and `loader/` (diagnostic fixtures for the binding runners) |
| `scripts/` | Validation, conformance orchestration, golden regeneration, thin per-binding runners |
| `docs/` | Hugo site; grid/rule/regridding pages are generated from the library files at build time |
| `archive/` | The pre-0.8.0 catalog, kept read-only until its migration completes |

## Quickstart

```sh
# Sibling checkout of EarthSciAST (or set ESS_ROOT):
git clone https://github.com/EarthSciML/EarthSciAST ../EarthSciAST

# Validate every library file (JSON schema + policy lint):
pip install jsonschema
python scripts/validate-library.py

# Run the cross-language conformance suite (needs the ESS binding toolchains):
./scripts/test-conformance.sh
```

## Status

The cross-binding conformance suite is green at **676 passed / 0 failed / 284
scope-skipped** across all five bindings, with `validate-library.py` reporting **0 findings
across 137 files**. The `ast` category is **byte-identical across all five bindings** (Julia,
Python, Rust, TypeScript, Go). The **numeric** categories (simulation, convergence,
regridding, reprojection) run on Julia (reference, `Tsit5`), Python (`LSODA`), and Rust; Go
and TypeScript are rewrite-only ports (`scope_excluded` from numeric).

Rust runs every numeric case its ODE solver can handle. The remaining nonlinear / high-order /
stiff cases are marked `blocked_upstream_bindings.rust` — for a specific reason established by
a **per-case tolerance sweep** (all three Rust solvers `Erk`/`Bdf`/`Sdirk` tested):
`earthsci-ast-rs`'s diffsol integrator has **no `dtmin`/max-step fail-fast guard**, so at
the tight tolerance the cross-binding gates require it drives `dt → 0` without terminating.
**33 cases** are blocked, in two honest classes — **22 hang at _every_ tolerance** (down to the
loosest `reltol 1e-6`: genuinely-stiff 2-D Laplacians, the `sqrt(0)`-singular Godunov
Hamiltonian, huge WENO ASTs), and **11 terminate only at a tolerance so loose** (`1e-6`–`1e-8`)
that their temporal error then misses the gate against Julia's tight reference (the `1e-9`
cross-binding-actuals floor for simulation, the `1e-4` errors-vs-golden for convergence). Two
low-order latlon zonal-advection *convergence* cases **do** run — at a pinned looser `reltol
1e-7`, loose enough not to hang yet tight enough that the spatially-dominated error still
matches the golden. Julia and Python verify every blocked case. The CI Rust job runs
`ast + numeric`; a blocked case skips cleanly (a `blocked-upstream` skip, never a silent
pass). All five ESS binding runners are registered in `scripts/test-conformance.sh` and CI.

Large AST goldens (≥ 64 KiB — the WENO / HJ-WENO reconstructions, whose fully-inlined
trees are 90–93 % redundant) are checked in as a committed **sha256 digest** rather than
the full bytes: every binding still reproduces the exact canonical AST and is gated on the
hash (`golden-digest`), so the repo carries no multi-MB derived blobs while keeping both
the byte-identity gate and the regression pin (AGENTS.md §5). (Shrinking these goldens
further — a canonical common-subexpression-elimination pass over the ~90 %-redundant lowered
trees — is an upstream ESM concern.)

The `ast` category is **byte-identical across Julia, Python, Rust, TypeScript,
and Go** — including the end-to-end consuming-model gate at N=64 and the two
newest spec mechanisms the library now leans on: §9.7.7 import
renaming/rebinding plus §9.6.1 `where` scoping (`two_cartesian_grids_coexist` —
one grid+rule imported twice, each instance rewritten to its own renamed axis
and shape) and §9.6.2 aggregate-mapped template expansion (`lcc_grid_roundtrip`
— reprojection templates inlined inside an `aggregate`).

Julia (reference), Python, and Rust run the numeric categories — MMS
simulation, convergence sweeps (error norms within rtol 1e-4 of the committed
Julia goldens, including the arbitrary-extent `heat_1d_zero_grad_nonunit`, the
parameterized fixed-value/fixed-flux BC drivers (`heat_1d_dirichlet`,
`heat_1d_neumann_flux`), the 2-D Laplacian heat driver (`heat_2d_neumann_flux`),
the centered-advection companion (`advection_1d_periodic_central`), the
fourth- and sixth-order periodic-heat drivers (`heat_1d_periodic_o4`,
`heat_1d_periodic_o6`) at observed order 4 and 6, the smoothly-stretched
non-uniform-mesh driver (`heat_1d_nonuniform`, conservative finite volume, L2
order 2 by supraconvergence), and both
lat-lon MMS drivers at observed order ~2), regridding, reprojection point
gates and the in-model LCC round-trip, and the ragged-MPAS divergence
simulation (Julia; div∘curl exact to ~3e-14). This wave adds the full
**first-derivative (gradient) family across all five grids** with second-derivative /
Laplacian companions, plus a **variable-coefficient / nonlinear Laplacian** `∇·(k∇u)` and
a **mixed `∂²/∂x∂y`** cross-derivative on cartesian — lifted to rank 3 on a non-uniform
vertical as the **`K_zz` boundary-layer mixing** operator with a **Robin surface-exchange**
(dry-deposition + emission) ground — the first **advection rule that matches the wind as an
operand**, `D(W·q, lev)`, taking a genuinely 3-D face-staggered velocity and picking its
donor from the *sign* of the local face velocity — the **metric spherical
(Laplace–Beltrami)** Laplacian on lat-lon, the **MPAS TRiSK edge-gradient and cell
Laplacian**, and a family of nonlinear high-order schemes: the **Godunov gradient-norm
Hamiltonian** (1-D and 2-D, exact on linear fields, entropy-fixed eikonal), fifth-order
**WENO-Z advection** and the **Jiang–Peng HJ-WENO** `|∇u|` (both observed order 5.00),
Colella–Woodward **PPM** conservative transport — reconstructed with the *full* CW84
limiter pair (eq. (1.8) monotonized slopes **and** the eq. (1.10) parabola limiter), so a
non-negative tracer with sharp jumps stays exactly within its initial range; using eq.
(1.10) alone, as this library formerly did, leaves the scheme unbounded and it manufactures
negative concentrations — and TVD **Lax–Friedrichs / minmod / superbee** limiters.
Sub-nominal observed orders (PPM's smooth-extremum clip — most visibly on the meridional
`cos⁴` case, whose smooth equatorial maximum the limiter flattens to first order — limiter
clipping, Lax–Friedrichs reducing to first-order upwind for a linear flux) are pinned as
*observed*, never forced to a design order. Regridding is now **end-to-end
declarative** — grid-spec → cell rings → geometry-derived broad-phase bin keys →
candidate-gated overlap → apply — with the gated path value-identical to the
dense path (tol 0.0) and conservation / partition-of-unity as the exact gates.
A **cartesian↔cartesian** remap constructs *both* grids in-library by importing
`cartesian_cell_rings.esm` twice under §9.7.7 prefixes (`src`/`tgt`), so neither
side is precomputed (`cartesian_rings_regrid_cart2cart_3x3_to_2x2`).

Scope gaps are recorded in the manifests, never shimmed: Go/TypeScript are
rewrite-only ports (`scope_excluded` from the numeric categories), and the
`blocked_upstream_bindings` entries name their precise upstream sites — notably
**Python spherical regridding stays gated on the optional `spherely` dependency**
(planar regrid runs on Julia + Python + Rust; spherical on Julia + Rust, and
activates for Python the moment the pinned wheel installs).

The repo carries exemplar content — uniform 1-D cartesian (arbitrary extent;
a centered second derivative with zero-gradient, parameterized Dirichlet, and
parameterized Neumann fixed-flux BCs; fourth- and sixth-order periodic Laplacians
that stay high-order to the wrap boundary; a centered first-derivative gradient; and
first-order upwind), uniform 2-D cartesian (the Laplacian assembled per-axis as
`D(D(u,x),x) + D(D(u,y),y)` with parameterized Neumann flux BCs), non-uniform 1-D
cartesian (a smoothly-stretched mesh from a single consumer-supplied edge array, with
the grid deriving cell centers and widths through its own `nonuniform_cell_center`/
`nonuniform_cell_width` templates over an `N+1`-node axis, feeding a conservative
finite-volume Laplacian with exact zero-flux walls, second-order in L2 by
supraconvergence), the lat-lon
production kit (coordinate/metric templates; periodic-lon and zero-gradient-lat
rules; global and regional recipes), its 3-D `latlon3d` extension for
GEOS-Chem-Classic-style tracer transport (importing the lat-lon horizontal
geometry and adding a hybrid non-uniform vertical through the Design-B
consumer-supplied edge-array pattern, with first-order-upwind, centered, and a
full rank-3 Colella–Woodward **PPM** family per direction — periodic-wrap zonal,
zero-gradient meridional, and a conservative **flux-form** vertical PPM whose
velocity vanishes at the model top/surface for exact machine-precision mass
conservation, offered in two variants: an *unlimited* reconstruction with verified
4th-order accuracy on the non-uniform mesh, and a **monotone** one (full CW84: eq
(1.6) edge on eq (1.8) monotonized slopes plus the eq (1.10) parabola limiter) that
is the production choice for chemistry-transport tracers — it holds a positive
tracer non-negative across sharp vertical gradients, where the unlimited scheme
measurably goes negative, at the cost of dropping to ~2.4 (L2) / ~1.9 (Linf) order;
both conserve mass to the bit — each verified by its own MMS convergence case, with
the upwind trio composing on one `[lon,lat,lev]` field into a full 3-D advection
driver; and, alongside transport, **vertical turbulent diffusion** `∂/∂z(K_zz ∂c/∂z)`
— the boundary-layer mixing operator — as a conservative non-uniform flux-difference
whose interior face fluxes telescope to machine-zero mass drift, offered with a
zero-flux lid plus either a zero-flux ground or a **Robin surface-exchange** ground
carrying dry deposition out and emission in. Because the diffusion rule matches the
*compound* `D(kz·D(c,lev),lev)` at priority 10 while advection matches a plain
`D(c,lev)`, the two fire on distinct terms and compose in one model — transport +
mixing + surface exchange is the minimum viable chemical-transport column. The
surface flux is closed by **resistance in series**, `(v_d·c₁ − E)·a/(a+v_d)` with
`a = 2K₁/dz₁`, and that is not a refinement but a correctness requirement: applying
the deposition velocity directly to the first-cell value — the obvious spelling — is
*inconsistent*, committing an O(dz) flux error that the divergence then divides by
`dz`, and it hides at coarse resolution (apparent order 1.98 at NLEV=32) before
unravelling (1.52 at 64, 1.16 at 128). The series form holds a clean 2.00, so the
surface exchange costs nothing in accuracy. Advection then takes the step the rest of the
transport stack is blocked on: **the wind becomes an operand of the rule**. The vertical PPM
above is already flux-form, but it reads its velocity as a *free name* — `w_edge`, a 1-D,
`lev`-only, time-static column profile — while every other advection rule (`upwind1_D_lon`,
the `central_D_*` family, zonal/meridional PPM) matches a *bare* `D(q,axis)` and has a
**scalar** velocity multiplied in from outside. That outside-multiplication is only correct
for a *constant* wind — `w·∂q/∂z` and `∂(w·q)/∂z` coincide only then — and silently stops
conserving tracer mass the moment the wind varies in space, as every real wind does. Worse,
a bare-`D` stencil has no velocity argument at all (`upwind1_D_lon_interior` is literally
`(f[i]−f[i−1])/dlon_deg`, a hard-coded backward difference), so it cannot switch donor when
the wind reverses. Matching the *compound* `D(W·q,lev)` instead binds the wind as a rule
parameter, which lets the operator take a real `[lon,lat,lev_nodes]` face-staggered field,
difference the flux `w·q` at faces so mass telescopes for *any* wind, and see `sign(w)` at
each face. The donor selection is written branch-free as `F = ½[w(q_L+q_R) − |w|(q_R−q_L)]`
— identical to `ifelse(w>0, w·q_L, w·q_R)` but with no boolean subexpression and no
removable singularity at `w=0`. Crucially, this compound match does **not** collide with the
diffusion rule even though `D(kz·D(u,lev),lev)` is *also* a `D` of a two-factor product: by
§9.6.1 a `where` shape constraint is satisfied only when the bound sub-AST is a **bare
variable reference**, and "a compound sub-AST fails it", so binding `q := D(u,lev)` fails the
constraint and the advection rule is filtered out at that node *before* priority selection
ever runs. The separation is a structural guarantee of the shape-constraint semantics, not a
priority race — verified by carrying both terms in one equation, where diffusion contributes
exactly zero on a constant field and advection exactly `−∂w/∂z` (agreement 2.2e-16), with the
column mass budget closing to 1.9e-16 relative, under one ulp. What this does *not* yet buy
is **free-stream preservation**: with a divergent wind a constant tracer does not stay
constant, because that requires consistency with continuity — the air-mass equation carried
alongside, `Δp·q` prognostic, with the *same* discrete face mass-fluxes driving both), and
the MPAS unstructured grid;
finite-difference/finite-volume rules; the conservative overlap regridder with
in-library cell-ring constructors; and Lambert conformal reprojection —
establishing the layering, testing, and docs patterns. The parameterized BCs
follow the reals-are-consumer-supplied contract: wall values/fluxes are free
names defaulting to 0, so the homogeneous case is the default and the Neumann
rule generalizes the zero-gradient one. A prototype **`duo`** icosahedral grid ships its
level-0 construction as pure closed-form AST (golden-ratio vertices normalized onto the
sphere via a nested aggregate) — a scoping result establishing that *resolution-parameterized*
subdivision needs two upstream ESM features (a `^`/pow metaparameter-expression op for the
`20·4^level` sizing and a build-time repeat/fold to iterate the refine pass, tracked as ESS
`ess-vnk`), with fixed levels otherwise shippable as MPAS-style const mesh data. The pre-0.8.0
catalog in `archive/` migrates rule-by-rule on top of these patterns.
