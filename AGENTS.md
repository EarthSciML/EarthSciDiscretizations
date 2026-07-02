# Agent and contributor guide

This repo is the ESM discretization standard library: pure data (`.esm` files) plus the
test harness and docs pipeline that pin its behavior. The sibling spec repo is
[EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization) (`$ESS_ROOT`,
default `../EarthSciSerialization` — resolved everywhere by `scripts/ess-locate.sh`).

## 1. AST first

All mathematics is expressed in the ESM Expression AST (esm-spec §4). If it fits on
paper, it fits in the AST. No helper functions hiding math, no registered functions for
anything writable in finite closed form (esm-spec §1.1).

## 2. Single pathway (ABSOLUTE)

This repo contains **zero evaluators, zero rule engines, zero numeric kernels**. Every
number in every golden was produced by an *official ESS binding runner* through the one
canonical pipeline: `.esm → parse → import/metaparameter resolution → §9.6.3 rewrite
fixpoint → official runner`.

- `scripts/runners/*` are thin argument-marshalling wrappers over ESS binding entry
  points. If a binding lacks a capability (e.g. a canonical-emit CLI), file the issue
  upstream in EarthSciSerialization — **never** patch a local shim, shadow evaluator, or
  per-rule special case here. The archive's retired shadow evaluators are the cautionary
  tale.
- `tools/` and `docs/generate_catalog.py` may *structurally* read `.esm` JSON (offsets,
  tags, metadata) and may call official ESS display/pretty-print APIs; they must never
  numerically evaluate library math.

## 3. Layering rules

- **Grid-first layout**: everything specific to a grid lives under `grids/<grid>/` —
  `grid.esm` (index sets + geometry, metaparameterized), `stencils/` (interior-only,
  match-less templates importing `grid.esm`), `rules/` (the complete auto-applied `match`
  rules: imported stencil + boundary-condition face regions in one `makearray`). The same
  scheme takes a different form on different grids; there is no grid-generic rule file.
- Boundary conditions live inside the rule body and nowhere else (esm-spec §9.6.8).
- Cross-grid entries (`regridding/`, `reprojection/`) and drivers (`problems/`) are
  top-level; problems import their grid + rules and carry inline tests (esm-spec §6.6.5).
- Everything is resolution-generic via metaparameters (esm-spec §9.7.6). Never bake a
  concrete size into a rule; bind sizes at import/subsystem edges (convergence wrappers)
  or defaults (the grid file).
- A rule author commits exactly: the library file(s), the `tests/conformance/ast/<rule>/`
  fixture + golden, a `tests/conformance/convergence/` case when an order claim is made,
  and nothing else.

## 4. Metadata and tags contract (lint-enforced)

`scripts/validate-library.py` enforces (codes L001–L008; fixtures in
`tests/invalid/lint/`):

- Exactly one kind tag: `esd:grid | esd:stencil | esd:rule | esd:regrid | esd:reproject |
  esd:problem` (L001), with required companions (L002):
  - `esd:grid` → `family:`
  - `esd:stencil` → `family: grid: op:`
  - `esd:rule` → `family: grid: op: order: bc: axes:`
  - `esd:regrid` → `method:` · `esd:reproject` → `crs:` · `esd:problem` → `grid:`
- Name agreement (L003): filename stem == `metadata.name` == the sole
  `expression_templates` key. Grid constructors are `grids/<name>/grid.esm` with
  `metadata.name == <name>`. Exception: `esd:regrid` / `esd:reproject` files factor
  one operation into several templates (weight/normalize/apply stages;
  forward/inverse pairs plus helper fragments), so they may declare multiple
  templates — every key must be the stem or prefixed `<stem>_`.
- `grid:<name>` resolves to `grids/<name>/grid.esm` (L004); stencils/rules live under
  that directory (L005).
- Rule `makearray` regions tile the axes named by `axes:` exactly, verified by folding
  metaparameter bounds at sampled sizes (L006). `axes:` lists the output dimensions in
  order, comma-separated (e.g. `axes:lon,lat`).
- Library entries are pure template-library files; problems carry models (L007). All
  files declare `esm: 0.8.0` (L008).
- Provenance goes in `metadata.references` (Fornberg, Snyder, Roache, …);
  human-readable derivations (incl. MMS solutions and BC-compatibility arguments) go in
  `metadata.description` / test descriptions.

## 5. Goldens

- Regeneration is **Julia-reference-only**, via `scripts/regenerate-goldens.sh
  [category …]` — it drives the canonical ESS pipeline and nothing else.
- A PR that changes goldens must say *why* (spec change, rule change, or bug — never
  "refreshed to green").
- Convergence goldens (`golden/errors.json`) hold the Julia-computed error norms;
  `scripts/check_convergence_order.py` asserts observed order from them without any
  binding installed. Docs read observed orders from the same files — nothing is ever
  recomputed for display.

## 6. Generated files

Docs content pages and plots are generated **only** in the CI docs job and are never
committed (`.gitignore` pins this). Conformance goldens and fixtures *are* committed.
Convergence per-resolution fixtures are thin hand-written wrappers binding
metaparameters — if you feel the need for a fixture generator, the metaparameter
mechanism is being misused.

## 7. Workflow

- Conventional commits: `type(scope): description` (scopes: `grids`, `rules`,
  `regridding`, `reprojection`, `problems`, `conformance`, `docs`, `ci`).
- Non-interactive shell flags for `cp`/`mv`/`rm`.
- CI: `.github/workflows/conformance.yml` (validate job always; per-binding matrix +
  compare as runners land), `docs.yml`, `nightly-upstream.yml` (drift watch vs ESS main).

## Phasing (bootstrap)

1. ✅ Skeleton + validation harness (this document, `validate-library.py`, lint fixtures)
2. ✅ Grids + stencils + rules + `ast/` conformance (Julia/Python runners)
3. ✅ Remaining AST runners (Rust/TypeScript/Go) — five-way byte-identical goldens
4. ✅ Problems + MMS + convergence (Julia reference + Python; Rust
   blocked-upstream per the manifests' `blocked_upstream_bindings` notes)
5. ✅ Regridding + reprojection (exact invariants + analytic point gates;
   per-pair regrid weights blocked-upstream per the case manifest)
6. ✅ Docs pipeline · 7. ✅ Full CI matrix (per-binding jobs + cross-binding
   compare in `.github/workflows/conformance.yml`)

The bootstrap is complete; the archive migration (`archive/` → `grids/` et al.)
proceeds rule-by-rule on top of these patterns. Each §9.7-capable binding runs
the categories in its scope (CONFORMANCE_SPEC.md §5.9); scope gaps are recorded
in the manifests (`scope_excluded`, `blocked_upstream_bindings`), never papered
over locally.
