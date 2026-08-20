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
its `metaparameters` (those *sizes*, with defaults), and geometry templates
that give the index-to-coordinate map (`x_coord`, `lon_coord`/`lat_coord`, the
`coslat` metric, …). The **real-valued extent is not baked in**: the origins
and spacings are *consumer-supplied free names* — `x0`/`dx` on
`cartesian_uniform_1d`, `lon0_deg`/`dlon_deg`/`lat0_deg`/`dlat_deg` (and
`R_sphere` for physical-units metrics) on `latlon` — that resolve in the
consuming model's scope at evaluation, the same keyed-factor contract MPAS uses
for `areaCell`/`dvEdge`. Only the integer *counts* are metaparameters; the
geometry itself is the consumer's, so one grid file serves any domain extent,
not merely any resolution (see [arbitrary extents](#arbitrary-extents-via-consumer-free-names)
below). On unstructured grids (MPAS) the whole geometry is mesh *data* — the
keyed factors `nEdgesOnCell`, `edgesOnCell`, `areaCell`, and (since the
regridding wave widened the contract) the vertex-geometry factors that let cell
polygon rings derive in-document.

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

## Arbitrary extents via consumer free names

Metaparameters give you every *resolution*; consumer-supplied free names give
you every *extent*. A grid file never fixes the physical domain — its geometry
templates and rule bodies reference bare names (`x0`, `dx`, `lon0_deg`, …) that
the consuming model must define as ordinary variables. On
`cartesian_uniform_1d`, for a domain `[a, b]`:

```json
"variables": {
  "x0": { "type": "parameter", "units": "m", "default": -1.5 },
  "dx": { "type": "unknown",   "units": "m" }
},
"equations": [
  { "lhs": "dx", "rhs": { "op": "/", "args": [4, "N"] } }
]
```

esm 1.0.0 declares exactly two variable types, `unknown` and `parameter`
(esm-spec §6.3.1). `dx` is an *observed* — but that is a **derived** class now,
not a declaration: a variable has no `expression` field, and what makes `dx`
observed is the bare-variable-LHS equation defining it. Ask the binding's
classification API (`observed_definitions`, `observed_unknowns`) when you need
to know; never read a declared type to answer it.

`dx = (b − a)/N` is written **dividing by the metaparameter name `N`**, which
§9.7.6 substitutes as an integer literal at load — so a convergence sweep that
rebinds `N` through the loader API keeps `dx` consistent automatically, with no
per-resolution files. `problems/heat_1d_zero_grad_nonunit.esm` proves it: the
*same* rule file that runs on the unit interval runs on `x ∈ [−1.5, 2.5]` at
observed order 2.00. A model that instead closes `N` at the import edge spells
the matching literal (`{op: /, args: [1, 8]}` for N = 8).

The `latlon` kit takes the same shape, region-generically: the **global** recipe
(`lon0_deg = 0`, `dlon_deg = 360/NLON`, `lat0_deg = −90`, `dlat_deg =
180/(NLAT − 1)`) closes the zonal circle so the periodic-lon rule applies; a
**regional** recipe (any other origin/spacing) does not close, so the
zero-gradient lon rule is used instead.

**Build-time scope caveat.** Free-name geometry templates (`x_coord`,
`lon_coord`) resolve only in *runtime* expression positions — rule bodies and
equation right-hand sides. They do **not** resolve in `ic` equations or §6.6.5
test `reference`s, whose build-time cellwise evaluation is scope-free in every
binding (model parameters/observeds are out of scope there). In those positions
spell the coordinate inline from literals and the metaparameter name — e.g.
`x_i = 0 + (i − 1/2)·(1/N)` — never `apply_expression_template x_coord`. See
[MMS and convergence](/guide/mms-and-convergence/).

## Two instances of one grid: rename and where

A rule fires by matching an operator (`D(D(f,x),x)`, `div(F)`, …). The `wrt`
literal and the operator alone are not enough to keep two grids apart when both
reuse an axis name, so every library rule also carries a §9.6.1 **`where`**
shape constraint — `{"f": {"shape": ["x"]}}`, `{"F": {"shape": ["edges"]}}`.
The rule then fires **only** on a bare field declared over *that* grid's index
sets, not on any unrelated derivative reusing the axis name.

Because the shape names an index set, it travels with §9.7.7 **import
renaming**. Import one grid+rule family twice under distinct prefixes and each
instance is scoped to its own mesh:

```json
"expression_template_imports": [
  { "ref": ".../rules/central_D2_zero_grad_bc.esm",
    "prefix": "meshA", "bindings": { "N": 8 },  "rebind": { "dx": "meshA_dx" } },
  { "ref": ".../rules/central_D2_zero_grad_bc.esm",
    "prefix": "meshB", "bindings": { "N": 16 }, "rebind": { "dx": "meshB_dx" } }
]
```

Renaming is transitive and simultaneous: the index set `x` arrives as
`meshA.x` / `meshB.x`, the match `wrt` follows it, `N` folds per edge into the
`makearray` regions, and each instance's `where` shape is rewritten to
`[meshA.x]` / `[meshB.x]`. The `rebind` retargets the free spacing name `dx` to
each mesh's own observed (`meshA_dx = 1/8`, `meshB_dx = 1/16`), so the two
fields evolve completely independently with no first-declared-wins collision and
no priority arbitration. `problems/two_cartesian_grids_coexist.esm` is the
worked driver — the composition that was inexpressible before renaming +
match-scoping landed. When you author a `where`-constrained rule, the consumer
must present a **bare declared shaped field**; a compound inline flux
(`div(u*h)`) must first be bound to a declared shaped observed.

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
coupling or model invokes with `apply_expression_template`.

Both are now **end-to-end declarative**. The conservative regridder ships its
whole pipeline in-library — cell-ring constructors (`cartesian_cell_rings.esm`
from a grid spec; MPAS `mpas_cell_rings_deg` from mesh vertex data), a broad
phase of geometry-derived integer `skolem` bin keys, and candidate-gated
overlap → normalize → apply templates whose gated form is value-identical to the
dense form under an admissible binning. Reprojection templates apply not only at
scalar points but **in-model over coordinate arrays**, by invoking the
forward/inverse templates inside an `aggregate` (§9.6.2 Option A expansion; the
`lcc_grid_roundtrip` case pins byte-identical lowering across all five
bindings). See the generated [Regridding](/regridding/) and
[Reprojection](/reprojection/) pages.
