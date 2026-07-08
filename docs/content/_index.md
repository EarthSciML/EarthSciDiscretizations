---
title: "EarthSciDiscretizations"
description: "The ESM discretization standard library: grids, finite-difference/finite-volume rewrite rules, regridding and reprojection expressions, with cross-language conformance goldens."
---

# EarthSciDiscretizations

The standard library of discretization rules for the
[EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization)
(ESS/ESM) format: grids, finite-difference/finite-volume rewrite rules,
regridding and reprojection expressions — together with the cross-language
conformance goldens, convergence suites, and MMS tests that pin their behavior
across every ESS binding.

This is a **pure data library**. Every entry is a `.esm` document; there is no
evaluator, rule engine, or numeric kernel in the repository. The ESM spec
(esm-spec.md §9.6.8) ships no discretization rules itself and delegates the
standard library here.

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

Boundary conditions live **inside** the rule (interior region = the stencil
aggregate, boundary-face regions = the BC); there is no separate
boundary-condition declaration anywhere (esm-spec §9.6.8). Rebinding the
metaparameters at the import edge is the whole convergence story — the same
files serve every resolution.

The physical *extent* is the consumer's too: a grid's origins and spacings are
supplied as free names (`x0`/`dx`, `lon0_deg`/`dlon_deg`/…) resolved in your
model's scope, so one grid file serves any domain, not just any resolution. And
because every rule carries a §9.6.1 `where` shape constraint that travels with
§9.7.7 import renaming, you can import one grid+rule family twice in a single
model — two meshes, each scoped to its own axis and spacing.

## Quickstart

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

Load the model through any official ESS binding: the import registers the rule,
the §9.6.3 rewrite fixpoint lowers `D(D(u,x),x)` into the rule's
`makearray` + `aggregate` stencil (boundary faces included), and the ordinary
simulation pathway takes it from there.

## Browse

- **[Grids](/grids/)** — one page per grid: index sets, geometry,
  stencils, and every rule published for it (match pattern, discretization,
  boundary treatment, observed convergence orders from the committed goldens).
- **[Regridding](/regridding/)** — conservative/interpolating remap expressions,
  end-to-end declarative: cell-ring constructors, an in-library broad phase, and
  candidate-gated overlap (gated == dense).
- **[Reprojection](/reprojection/)** — coordinate-transform template fragments,
  usable in-model over coordinate arrays via aggregates.
- **[Guide](/guide/)** — the layering, authoring a rule, MMS + convergence
  testing, and the conformance suite.
- Facets: [families](/families/), [operators](/ops/),
  [boundary conditions](/bcs/), [tags](/tags/).

The grid, regridding, and reprojection pages are **generated from the library
files at build time** by `docs/generate_catalog.py` (they are never committed);
the source of truth is always the `.esm` file linked at the top of each page.
