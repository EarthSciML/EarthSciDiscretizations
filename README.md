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

The conformance infrastructure is complete. All five ESS binding runners are
registered in `scripts/test-conformance.sh` and CI: the `ast` category is
**byte-identical across Julia, Python, Rust, TypeScript, and Go** (including the
end-to-end consuming-model gate at N=64); Julia (reference) and Python run every
category — MMS simulation, convergence sweeps (error norms within rtol 1e-4 of
the committed Julia goldens), regridding invariants, and reprojection point
gates. Scope gaps are recorded in the manifests, never shimmed: Go/TypeScript
are rewrite-only ports (`scope_excluded`), and the Rust §6.6.5 simulation
pathway is `blocked_upstream_bindings` pending its `run_pde_tests` port.

The repo carries exemplar content (uniform 1-D cartesian, lat-lon, and MPAS
grids; centered/upwind rules; a conservative regridder; Lambert conformal
reprojection) establishing the layering, testing, and docs patterns. The
pre-0.8.0 catalog in `archive/` migrates rule-by-rule on top of these patterns.
