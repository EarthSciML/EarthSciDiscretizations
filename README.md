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

The conformance suite is green at **192 passed / 0 failed / 46 scope-skipped**,
with `validate-library.py` reporting 0 findings. All five ESS binding runners
are registered in `scripts/test-conformance.sh` and CI.

The `ast` category is **byte-identical across Julia, Python, Rust, TypeScript,
and Go** — including the end-to-end consuming-model gate at N=64 and the two
newest spec mechanisms the library now leans on: §9.7.7 import
renaming/rebinding plus §9.6.1 `where` scoping (`two_cartesian_grids_coexist` —
one grid+rule imported twice, each instance rewritten to its own renamed axis
and shape) and §9.6.2 aggregate-mapped template expansion (`lcc_grid_roundtrip`
— reprojection templates inlined inside an `aggregate`).

Julia (reference), Python, and Rust run the numeric categories — MMS
simulation, convergence sweeps (error norms within rtol 1e-4 of the committed
Julia goldens, including the arbitrary-extent `heat_1d_zero_grad_nonunit` and
both lat-lon MMS drivers at observed order ~2), regridding, reprojection point
gates and the in-model LCC round-trip, and the ragged-MPAS divergence
simulation (Julia; div∘curl exact to ~3e-14). Regridding is now **end-to-end
declarative** — grid-spec → cell rings → geometry-derived broad-phase bin keys →
candidate-gated overlap → apply — with the gated path value-identical to the
dense path (tol 0.0) and conservation / partition-of-unity as the exact gates.

Scope gaps are recorded in the manifests, never shimmed: Go/TypeScript are
rewrite-only ports (`scope_excluded` from the numeric categories), and the
`blocked_upstream_bindings` entries name their precise upstream sites — notably
**Python spherical regridding stays gated on the optional `spherely` dependency**
(planar regrid runs on Julia + Python + Rust; spherical on Julia + Rust, and
activates for Python the moment the pinned wheel installs).

The repo carries exemplar content — uniform 1-D cartesian (arbitrary extent),
the lat-lon production kit (coordinate/metric templates; periodic-lon and
zero-gradient-lat rules; global and regional recipes), and the MPAS unstructured
grid; centered/upwind/finite-volume rules; the conservative overlap regridder
with in-library cell-ring constructors; and Lambert conformal reprojection —
establishing the layering, testing, and docs patterns. The pre-0.8.0 catalog in
`archive/` migrates rule-by-rule on top of these patterns.
