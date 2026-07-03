# EarthSciDiscretizations

The standard library of discretization rules for the
[EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization) (ESS/ESM)
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
# Sibling checkout of EarthSciSerialization (or set ESS_ROOT):
git clone https://github.com/EarthSciML/EarthSciSerialization ../EarthSciSerialization

# Validate every library file (JSON schema + policy lint):
pip install jsonschema
python scripts/validate-library.py

# Run the cross-language conformance suite (needs the ESS binding toolchains):
./scripts/test-conformance.sh
```

## Status

The **Julia reference runner passes all 140 cases** (45 AST · 46 simulation · 40
convergence · 8 regridding · 1 reprojection) and `validate-library.py` reports **0
findings across 137 files**. The `ast` category is **byte-identical across all five
bindings** (Julia, Python, Rust, TypeScript, Go — 405 passed / 0 failed), and the
numeric categories are cross-verified **Julia ↔ Python** (298 passed / 0 failed / 61
scope-skipped — the `spherely`-gated spherical regrids plus the Go/TypeScript numeric
scope). In CI the **AST** gate covers all five bindings, but the **numeric** categories
run on **Julia and Python only**. The Rust CI job is AST-only for a concrete upstream
reason: `earthsci-toolkit-rs` offers only the explicit **`Erk`** integrator, which **hangs
on genuinely-two-axis stiff diffusion** (confirmed on `anisotropic_diffusion_2d_periodic`,
whose full `D2x + Dxy + D2y` tensor is stiff in both directions at once) — an explicit RK
method cannot step a stiff 2-D Laplacian. Rust completes the 1-D and hyperbolic numeric
cases (the per-axis `heat_2d_*` drivers diffuse only one axis and behave like the 1-D
heat cases), but a full Rust numeric sweep hangs, so it is not wired into CI (tracked as
an upstream stiff-solver gap). Julia (`Tsit5`) and Python (`LSODA`) verify every numeric
case; Go/TypeScript are rewrite-only ports. All five ESS binding runners are registered in
`scripts/test-conformance.sh` and CI.

Large AST goldens (≥ 64 KiB — the WENO / HJ-WENO reconstructions, whose fully-inlined
trees are 90–93 % redundant) are checked in as a committed **sha256 digest** rather than
the full bytes: every binding still reproduces the exact canonical AST and is gated on the
hash (`golden-digest`), so the repo carries no multi-MB derived blobs while keeping both
the byte-identity gate and the regression pin (AGENTS.md §5). A companion upstream RFC
proposes an ESM `let`/shared-binding node + canonical CSE pass to shrink such goldens ~10×
at the source (`EarthSciSerialization/docs/content/rfcs/shared-subexpression-binding-cse.md`).

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
a **mixed `∂²/∂x∂y`** cross-derivative on cartesian, the **metric spherical
(Laplace–Beltrami)** Laplacian on lat-lon, the **MPAS TRiSK edge-gradient and cell
Laplacian**, and a family of nonlinear high-order schemes: the **Godunov gradient-norm
Hamiltonian** (1-D and 2-D, exact on linear fields, entropy-fixed eikonal), fifth-order
**WENO-Z advection** and the **Jiang–Peng HJ-WENO** `|∇u|` (both observed order 5.00),
Colella–Woodward **PPM** conservative transport, and TVD **Lax–Friedrichs / minmod /
superbee** limiters. Sub-nominal observed orders (PPM's smooth-extremum clip, limiter
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
rules; global and regional recipes), and the MPAS unstructured grid;
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
