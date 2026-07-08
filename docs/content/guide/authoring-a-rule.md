---
title: "Authoring a rule"
description: "Walkthrough: grid file → stencil → BC rule → lint → fixtures → goldens, per the AGENTS.md contract."
weight: 2
---

A rule author commits exactly: the library file(s), the
`tests/conformance/ast/<rule>/` fixture + golden, a
`tests/conformance/convergence/` case when an order claim is made — and nothing
else (AGENTS.md §3). Docs pages are generated; never commit them.

## 0. Ground rules

- **AST first.** All mathematics is the ESM Expression AST (esm-spec §4). No
  helper functions hiding math, no registered functions for anything writable
  in finite closed form.
- **Single pathway.** This repo contains zero evaluators. Every golden number
  is produced by an official ESS binding runner through
  `.esm → parse → import/metaparameter resolution → §9.6.3 rewrite fixpoint →
  official runner`. If a binding lacks a capability, file the issue upstream in
  EarthSciAST — never patch a local shim.
- Everything specific to a grid lives under `grids/<grid>/`. There is no
  grid-generic rule file: the same scheme takes a different form on different
  grids.

## 1. The grid file (reuse it if it exists)

`grids/<grid>/grid.esm` declares metaparameters, index sets, and geometry
templates, and carries the `esd:grid` kind tag plus `family:`. The filename is
`grid.esm` and `metadata.name` equals the directory name (lint L003).

## 2. The stencil

`grids/<grid>/stencils/<name>.esm`: a match-less template importing
`../grid.esm`, whose body is the interior `aggregate` only. Tags:
`esd:stencil` + `family:`, `grid:`, `op:` (L002). Put the provenance (Fornberg,
LeVeque, Ringler, …) in `metadata.references` and the derivation in
`metadata.description`.

```json
"central_D2_interior": {
  "params": ["f"],
  "body": {
    "op": "aggregate", "output_idx": ["i"], "args": ["f"],
    "ranges": { "i": [2, { "op": "-", "args": ["N", 1] }] },
    "expr": { "...": "(f[i+1] − 2·f[i] + f[i−1]) / dx²" }
  }
}
```

## 3. The rule

`grids/<grid>/rules/<name>.esm`: imports the stencil (or the grid directly for
one-region schemes), declares the single template whose `match` fires on the
target operator and whose body is one `makearray`:

- **regions** tile the axes named by the `axes:` tag exactly — lint L006
  verifies this by folding metaparameter bounds at sampled sizes (one even,
  one odd — parity mistakes in region arithmetic are the classic bug).
- The **interior region** invokes the stencil via
  `apply_expression_template`.
- The **face regions** are the boundary condition: one-sided differences for
  zero-gradient, a periodic wrap, a fixed value for Dirichlet, an `index` into
  another variable for a shared seam. Later regions overwrite earlier
  (esm-spec §4.3.2).
- A compound scheme (e.g. matching `D(D(f,x),x)` whole) sets `priority` so it
  out-ranks plain-D rules (§9.6.3).
- A **`where` shape constraint** (§9.6.1) scopes the match to *this* grid:
  `"where": { "f": { "shape": ["x"] } }` (or `{"F": {"shape": ["edges"]}}` for
  the MPAS div rule). Every library rule carries one — it stops an axis-name
  collision from firing the rule on an unrelated derivative, and it is what lets
  the rule be imported twice (renamed) without a first-declared-wins tie-break:
  the shape names an index set, so §9.7.7 renaming rewrites it to `[meshA.x]` /
  `[meshB.x]` along with the `wrt` literal. See
  [the layering](/guide/layering/#two-instances-of-one-grid-rename-and-where).

  **Authoring convention:** a `where` shape requires a *bare declared shaped
  field*. Write the constraint over a plain parameter (`f`, `F`), and document
  that a consumer differentiating a compound inline expression (`div(u*h)`,
  `D(D(u*v,x),x)`) must first bind it to a declared shaped observed — the
  constraint will not match a compound inline argument.

Tags: `esd:rule` + `family:`, `grid:`, `op:`, `order:`, `bc:`, `axes:` (L002);
`axes:` lists the output dimensions in order, comma-separated. The filename
stem, `metadata.name`, and the sole `expression_templates` key must agree
(L003), and the file must live under the grid directory its `grid:` tag names
(L004/L005). All files declare `esm: 0.8.0` (L008).

## 4. Lint

```sh
python scripts/validate-library.py            # schema + policy lint L001–L008
python scripts/validate-library.py --invalid  # lint-fixture expectations
```

This needs a sibling `EarthSciAST` checkout (or `ESS_ROOT`) for the
JSON schema, and `pip install jsonschema`.

## 5. Fixtures

- **`tests/conformance/ast/<rule>/`** — the post-lowering AST fixture: a tiny
  consuming model binding small sizes, plus the canonical-emit golden. Every
  binding must reproduce it byte-identically.
- **`tests/conformance/convergence/<case>/manifest.json`** — required whenever
  the rule claims an `order:`. It names the problem, the rule, the resolutions
  (loader-API bindings — no per-resolution files), norms, `expected_order`,
  and `order_tolerance`, and starts life with `"status": "pending-golden"`.
- A driving problem under `problems/` with inline §6.6.5 MMS tests, if one
  does not already exist (see [MMS and convergence](/guide/mms-and-convergence/)).

## 6. Goldens

Golden regeneration is Julia-reference-only:

```sh
./scripts/regenerate-goldens.sh convergence
```

It drives the canonical ESS pipeline and nothing else. A PR that changes
goldens must say *why* (spec change, rule change, or bug — never "refreshed to
green"). Once `golden/errors.json` lands, `scripts/check_convergence_order.py`
asserts the observed order in CI without any binding installed, and the docs
build renders the observed-order table and convergence plot from the same file.

## 7. Done

Open the PR with a conventional-commit title (`feat(rules): …`). CI runs the
validate job (schema + lint + convergence-order check + spellcheck), the
per-binding conformance matrix as runners land, and the docs build — your rule
appears on its grid's page automatically.
