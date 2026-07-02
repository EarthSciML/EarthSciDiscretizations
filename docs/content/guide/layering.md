---
title: "The layering"
description: "The compositional design: grid → stencil → rule → problem, held together by template imports and load-time metaparameters."
weight: 1
---

Every library entry is a valid ESM **template-library file** (esm-spec §9.7.1):
a document whose payload is `expression_templates`, plus optionally
`index_sets`, `metaparameters`, and `expression_template_imports`. Files
compose by reference — each layer imports the one below it — and nothing in the
chain hard-codes a resolution.

## The four layers

```
grids/<grid>/grid.esm             the grid constructor
  └─ stencils/<name>.esm          interior operator
       └─ rules/<name>.esm        stencil + boundary conditions
            └─ problems/<name>.esm   runnable driver (or your model)
```

**`grid.esm`** owns the axes and the geometry. It declares the grid's
`index_sets` (e.g. `x` of size `N`, or `lon`/`lat` of sizes `NLON`/`NLAT`),
its `metaparameters` (those sizes, with defaults), and zero-parameter geometry
templates — `dx = 1/N`, cell-center coordinates, spacings in degrees. Geometry
is ordinary expressions, resolution-generic through the metaparameters; on
unstructured grids (MPAS) geometry is mesh *data* instead, and the grid file
documents the keyed-factor contract a consumer must satisfy.

**`stencils/<name>.esm`** is the interior operator only: a **match-less** named
template importing `grid.esm`, e.g. the 3-point centered second difference over
`i ∈ [2, N−1]`. It cannot fire on its own — it exists to be wrapped.

**`rules/<name>.esm`** is the complete, auto-applied rewrite rule: a template
with a `match` pattern on a spatial `D` (or on `div`, `laplacian`, … — the open
namespace of esm-spec §4.2) whose body is a single `makearray`. The interior
region invokes the imported stencil; the boundary-face regions encode the
boundary condition. **Boundary conditions live inside the rule and nowhere
else** (esm-spec §9.6.8): a discretized derivative over a finite domain is
inseparable from its boundary treatment, so there is no separate BC declaration
anywhere in the format. Choosing a scheme = choosing a rule
(`central_D2_zero_grad_bc`, `upwind1_D_periodic`, …). Compound schemes out-rank
plain per-derivative rules via `priority` (esm-spec §9.6.3), so a
second-derivative rule at priority 10 fires on `D(D(f,x),x)` before any plain-D
rule lowers the inner derivative.

**`problems/<name>.esm`** is a runnable driver: a model importing its grid's
rule, with equations, initial conditions written over the grid's geometry
templates, and inline §6.6.5 MMS tests. Your own model sits at exactly this
layer — one import, one `D`.

## Metaparameters: one file, every resolution

Metaparameters are document-scoped, load-time integers (esm-spec §9.7.6). They
are admissible wherever the schema wants a structural integer — index-set
sizes, `aggregate` range bounds, `makearray` region bounds — and those sites
fold to concrete integers at load. In ordinary expression positions a
metaparameter name is substituted as an integer literal (`360/NLON` stays an
AST division).

Bindings flow down the reference DAG; open metaparameters flow up:

1. **Import edge** — `expression_template_imports[k].bindings` closes the
   imported document's metaparameters (`{ "N": 128 }`).
2. **Re-export** — metaparameters left unbound at an edge are inherited by the
   importer, so binding once at the top of a chain reaches the grid file at
   the bottom.
3. **Subsystem edge** — a §4.7 `{"ref": …}` may carry the same `bindings`
   (e.g. a convergence wrapper instantiating a problem at a given size).
4. **Loader API** — hosts may bind the root document's open metaparameters at
   load. This is how convergence sweeps run one problem file at N = 16, 32,
   64, 128 with **no per-resolution fixture files**.
5. **Defaults, last** — a still-open metaparameter takes its declared
   `default`.

Never bake a concrete size into a rule. If you feel the need for a fixture
generator, the metaparameter mechanism is being misused (AGENTS.md §6).

## What resolves when

Per document, innermost-first (esm-spec §9.7.6): resolve imports (instantiating
metaparameters at each edge) → merge `index_sets` into the importing document
→ fold this document's metaparameters → inline `apply_expression_template`
body references (§9.7.3) → run the §9.6.3 rewrite fixpoint on fully-concrete
trees. By the time any `match` rule is considered, every rule body is a closed
Expression AST. Two bindings expanding the same file must produce
byte-identical post-lowering ASTs — that determinism contract is what the
`ast/` conformance category pins across all five bindings.

## Cross-grid entries

`regridding/` and `reprojection/` are top-level because they connect grids
rather than belonging to one. They follow the same template-library form:
match-less templates (overlap-area matrices, projection formulas) that a
coupling or model invokes with `apply_expression_template`. See the generated
[Regridding](/regridding/) and [Reprojection](/reprojection/) pages.
