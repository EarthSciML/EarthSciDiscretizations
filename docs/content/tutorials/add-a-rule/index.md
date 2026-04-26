---
title: "Add a new discretization rule"
slug: "add-a-rule"
description: "End-to-end contributor walkthrough: declare a rule JSON, ship canonical + convergence fixtures, register with the walker, and render the catalog page. Worked example: centered_2nd_uniform on a Cartesian axis."
tags: ["tutorial", "contributor", "finite-difference"]
---

This tutorial walks through every step required to land a new discretization
rule in the ESD catalog, end-to-end. We use
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}) — the
two-point centered finite difference on a uniform Cartesian axis — as the
running example, because it is the smallest rule in the repo that exercises
the full pipeline (rule JSON → canonical fixture → MMS convergence sweep →
walker registration → doc page).

The audience is a contributor who is comfortable with Julia/MTK but new to
this repo. By the end of the tutorial you will know:

1. Where each artifact lives and what owns its schema.
2. How the three CI layers (A canonical-form round-trip, B MMS convergence,
   B′ monotonicity) are wired up.
3. How to land a new rule with the staged workflow we use in practice
   (`applicable: false` first, flip to `true` only after Layer B passes
   locally).

Throughout, file paths are clickable into the canonical example committed at
HEAD — open each one as you read.

## What a "rule" is in this repo

A **discretization rule** maps a continuous PDE operator (e.g.
`grad(u, dim=x)`) onto a discrete stencil over neighbors of a target cell,
plus a coefficient AST evaluated against per-grid bindings. Rule files live
under [`discretizations/<family>/`]({{< param repoURL >}}/tree/main/discretizations)
as JSON, and CI validates them against the EarthSciSerialization §7 schema.

The authoring policy is **AST-first**: every coefficient — including limiters
and reconstructions — is expressed directly in the
`ExpressionNode` AST that ESS already understands (`+ - * / ^`,
`max min ifelse abs sign`, all trig/exp/log, comparisons, logical ops).
ESD does **not** carry a shadow Julia evaluator; the only supported entry
point is
[`EarthSciDiscretizations.eval_coeff`]({{< param repoURL >}}/blob/main/src/rule_eval.jl),
a thin passthrough to `EarthSciSerialization.evaluate`. See
[`AGENTS.md`]({{< param repoURL >}}/blob/main/AGENTS.md) for the full
authoring policy.

## Prerequisites

Resolve the `EarthSciSerialization` dependency (it is not in the General
registry yet) and verify the package loads:

```bash
scripts/setup_polecat_env.sh
julia --project=. -e 'using EarthSciDiscretizations'
```

The script prefers a workspace checkout under
`../EarthSciSerialization/` and falls back to a `Pkg.add(url=...)` from
GitHub. It is idempotent — re-run any time.

## The running example

The rule we will reproduce is the two-point centered second-order finite
difference for ∂u/∂x on a uniform Cartesian axis:

$$\left(\frac{\partial u}{\partial x}\right)_i \;\approx\; \frac{u_{i+1} - u_{i-1}}{2\,\Delta x}.$$

Stencil: two neighbors at offset ±1 with symmetric coefficients
±1 / (2 dx). Grid family: `cartesian`. Combine: `+`. Accuracy:
O(dx²).

When you adapt this tutorial to a new rule, swap each artifact below for
the equivalent in your scheme — the layout and ordering of steps is the
same regardless of family.

## Step 1 — Declare the rule JSON

**Path:** `discretizations/<family>/<rule_name>.json`
**Reference:** [`discretizations/finite_difference/centered_2nd_uniform.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform.json)

The rule JSON has four pieces:

```json
{
  "discretizations": {
    "centered_2nd_uniform": {
      "applies_to":  { "op": "grad", "args": ["$u"], "dim": "$x" },
      "grid_family": "cartesian",
      "combine":     "+",
      "accuracy":    "O(dx^2)",
      "stencil":     [ /* one entry per neighbor */ ]
    }
  }
}
```

- **`applies_to`** — the *applicability declaration*. This is the AST pattern
  that the ESS rewriter matches against the equation tree. `$u` and `$x`
  are pattern metavariables that bind to the actual field name and axis at
  rewrite time. If your rule applies to multiple operators, write
  one rule file per pattern (the catalog favors small, single-purpose
  rules).
- **`grid_family`** — exactly one of `cartesian`, `vertical`, `latlon`,
  `cubed_sphere`, `mpas`, `duo`. The selector `kind` inside `stencil`
  entries is per-family; see
  [`discretizations/SELECTOR_KINDS.md`]({{< param repoURL >}}/blob/main/discretizations/SELECTOR_KINDS.md)
  decision #1 for why we did *not* go structural.
- **`combine`** — how the per-neighbor contributions reduce to the
  discrete operator. Almost always `+`.
- **`stencil[].selector`** — `{ kind, axis, offset }` for structured 1D
  rules. For MPAS use `kind: "indirect"` or `kind: "reduction"` (see
  `SELECTOR_KINDS.md` decisions #6–#8). For 2D in-panel rules
  (`covariant_laplacian_cubed_sphere`) each entry carries a `selectors`
  array.
- **`stencil[].coeff`** — an `ExpressionNode` AST. Available bindings
  depend on the grid family: cartesian binds `dx`/`h`, latlon binds
  `R`, `dlon`, `dlat`, `cos_lat`, etc. The full per-family symbol set is
  documented in `SELECTOR_KINDS.md`.

For our running example the two stencil entries are:

```json
{ "selector": { "kind": "cartesian", "axis": "$x", "offset": -1 },
  "coeff":    { "op": "/", "args": [-1, { "op": "*", "args": [2, "dx"] }] } },
{ "selector": { "kind": "cartesian", "axis": "$x", "offset":  1 },
  "coeff":    { "op": "/", "args": [ 1, { "op": "*", "args": [2, "dx"] }] } }
```

> **Naming convention.** When the same scheme exists for multiple grid
> families, append a per-family suffix: `centered_2nd_uniform.json` (cartesian),
> `centered_2nd_uniform_vertical.json`, `centered_2nd_uniform_latlon.json`.
> One file per (scheme, family) pair — see
> [`discretizations/finite_difference/README.md`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/README.md).

## Step 2 — Add the canonical fixture (Layer A)

**Path:** `discretizations/<family>/<rule_name>/fixtures/canonical/{input.esm, expected.esm}`
**Reference:** [`discretizations/finite_difference/centered_2nd_uniform/fixtures/canonical/`]({{< param repoURL >}}/tree/main/discretizations/finite_difference/centered_2nd_uniform/fixtures/canonical)

Layer A drives the rule end-to-end through
`EarthSciSerialization.discretize` and byte-compares the canonical-form
output against `expected.esm`. This is the conformance signature that
guarantees rewriter behavior is stable across bindings.

`input.esm` is a complete `.esm` document — grids + rules + a model with one
equation that contains the operator your rule matches. The minimal shape:

```json
{
  "esm": "0.2.0",
  "metadata": { "name": "<rule_name>_canonical" },
  "grids": { "gx": { "family": "cartesian", "dimensions": [ ... ] } },
  "rules": [ { "name": "<rule_name>", "pattern": ..., "replacement": ... } ],
  "models": {
    "M": {
      "grid": "gx",
      "variables": { "u": { "type": "state", ... } },
      "equations": [
        { "lhs": { "op": "D", "args": ["u"], "wrt": "t" },
          "rhs": { "op": "grad", "args": ["u"], "dim": "i" } }
      ]
    }
  }
}
```

`expected.esm` is the **canonical-form** rendering of `discretize(input)`:
sorted keys, minified, `format_canonical_float` for floats. The walker
compares byte-for-byte (trailing newlines tolerated), so do not hand-edit
this file beyond what the canonical emitter produces. To regenerate
locally:

```julia
julia --project=. -e '
  using JSON, EarthSciSerialization
  doc = JSON.parse(read("discretizations/finite_difference/<rule_name>/fixtures/canonical/input.esm", String))
  out = EarthSciSerialization.discretize(doc)
  println(EarthSciSerialization.canonical_doc_json(out))
' > discretizations/finite_difference/<rule_name>/fixtures/canonical/expected.esm
```

(The exact helper name lives in
[`test/walk_esd_tests.jl`]({{< param repoURL >}}/blob/main/test/walk_esd_tests.jl)
under `apply_rule_and_diff`/`canonical_doc_json` — copy that emission path
verbatim so the byte signature matches the walker.)

> **When to use `rewrite/` instead of `canonical/`.** Index-rewrite rules
> like `periodic_bc` operate on an expression rather than a whole document,
> so they ship a sibling `fixtures/rewrite/{input,expected}.esm` pair
> instead. The walker AND-combines both variants when present; see
> `apply_rewrite_and_diff` in `test/walk_esd_tests.jl`.

## Step 3 — Add the convergence fixture (Layer B), `applicable: false` first

**Path:** `discretizations/<family>/<rule_name>/fixtures/convergence/{input.esm, expected.esm}`
**Reference (working sweep):** [`discretizations/finite_difference/centered_2nd_uniform/fixtures/convergence/`]({{< param repoURL >}}/tree/main/discretizations/finite_difference/centered_2nd_uniform/fixtures/convergence)
**Reference (`applicable: false` form):** [`discretizations/finite_volume/flux_limiter_minmod/fixtures/convergence/`]({{< param repoURL >}}/tree/main/discretizations/finite_volume/flux_limiter_minmod/fixtures/convergence)

Layer B dispatches to ESS's `verify_mms_convergence`, which evaluates every
stencil coefficient through the AST evaluator and runs a manufactured-solution
sweep across a sequence of grids. The fixture pair is small but specific:

`input.esm`:

```json
{
  "rule": "<rule_name>",
  "manufactured_solution": "sin(2*pi*x) on [0,1] periodic; derivative 2*pi*cos(2*pi*x)",
  "sampling": "cell_center",
  "grids": [ { "n": 16 }, { "n": 32 }, { "n": 64 }, { "n": 128 } ]
}
```

`expected.esm`:

```json
{
  "rule": "<rule_name>",
  "metric": "Linf",
  "expected_min_order": 1.9,
  "notes": "Theoretically O(dx^2); 1.9 tolerates pre-asymptotic drift on the 16->32->64->128 sequence."
}
```

**Stage the fixture as not-applicable initially.** Until you have actually
run the sweep locally and seen the slope match theory, both files should
declare `applicable: false` so the walker reports a clean SKIP rather than
a possibly-wrong PASS:

```json
{ "rule": "<rule_name>", "applicable": false,
  "skip_reason": "Pending local Layer B verification — convergence fixture is staged but not yet run." }
```

This is the same idiom that
[`flux_limiter_minmod`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_minmod/fixtures/convergence/input.esm)
uses permanently (limiters do not have a convergence acceptance) and that
[`centered_2nd_uniform_latlon` used during landing]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_latlon/fixtures/convergence/expected.esm)
before its harness extension landed. The walker reports SKIP with reason
`"fixture-declared not applicable"` for these.

## Step 4 — Run Layer A to verify the canonical round-trip

```bash
julia --project=. test/runtests.jl
```

The harness uses `TestItemRunner` (see
[`test/runtests.jl`]({{< param repoURL >}}/blob/main/test/runtests.jl)).
To filter to just the walker test while iterating, use:

```bash
julia --project=. -e '
  using TestItemRunner
  @run_package_tests filter=ti->occursin("walker", ti.name)'
```

The walker test (in
[`test/test_esd_walker.jl`]({{< param repoURL >}}/blob/main/test/test_esd_walker.jl))
emits a JUnit XML at `test/junit-esd.xml` that lists Layer A / B / B′ / C
outcome per rule with a one-line reason. Open it after a run to confirm
your new rule shows up with `layer_a == LAYER_PASS`.

> **What "PASS" means here.** Layer A passes only when
> `EarthSciSerialization.discretize(input.esm)`'s canonical output matches
> `expected.esm` byte-for-byte. A whitespace-only difference is still a
> failure — regenerate `expected.esm` from the canonical emitter rather
> than hand-editing.

## Step 5 — Verify Layer B locally and flip `applicable: true`

Run only the walker (above) and check `test/junit-esd.xml` for your rule's
Layer B outcome. The reason field on PASS reads
`"min order ≈ <observed_order>"`; on FAIL it shows the per-grid residuals.

Once the sweep passes locally with comfortable margin, replace
`applicable: false` with the real metric and threshold:

```json
{
  "rule": "<rule_name>",
  "metric": "Linf",
  "expected_min_order": 1.9,
  "notes": "..."
}
```

Pick `expected_min_order` slightly below the theoretical order (e.g. 1.9
for an O(dx²) scheme) so minor pre-asymptotic drift on the 16 → 32 → 64 →
128 sequence does not flap the build. Document the slack in `notes`.

## Step 6 — Register the rule with the walker harness

**File:** [`test/test_esd_walker.jl`]({{< param repoURL >}}/blob/main/test/test_esd_walker.jl)

The walker discovers rules automatically via `load_rules`, but the test
asserts an explicit *expected* set so accidental deletions surface as
test failures rather than silent no-ops. Two edits:

1. Add your rule name to the `seeded` tuple near the top of the test (the
   "every rule we expect to find on disk" assertion).
2. Add the `(family, rule_name)` pair to one of:
   - `pass_layer_b` — Layer B is applicable and runs to PASS
   - `not_applicable_layer_b` — Layer B is permanently `applicable: false`
     (limiters, rules pending ESS harness extensions, BCs)

If your rule ships a `monotonicity/` fixture (Layer B′), also add it to
`pass_layer_limiter`. Layer C (integration benchmarks) is gated on
`ESD_RUN_INTEGRATION=1` and skipped by default.

## Step 7 — Catalog row + doc page

1. **Catalog row.** Add a row for your rule in
   [`docs/rule-catalog.md`]({{< param repoURL >}}/blob/main/docs/rule-catalog.md).
   Match the column shape used by the surrounding rows in your section
   (Section A for MOL parity, B for CAM5 FV core, etc.). Mark the
   `audit_status` as `looks_ok` for green-field rules, `needs_review` if
   you ported numerics from pre-Gas-Town `src/`.
2. **Hugo page.** Add `docs/content/rules/<rule_name>.md`. Use
   [`docs/content/rules/centered_2nd_uniform.md`]({{< param repoURL >}}/blob/main/docs/content/rules/centered_2nd_uniform.md)
   as a template — the frontmatter taxonomies (`families`, `grid_families`,
   `rule_kinds`, `tags`) drive the faceted navigation. The body sections
   we use are: **Stencil** (with the auto-generated PNG), **Coefficients**
   (table), **Discrete operator** (KaTeX block), **Convergence** (with the
   auto-generated PNG, or a *pending* callout if Layer B is not applicable).
3. **Plot artifacts.** Register your rule in
   [`tools/render_doc_plots.py`]({{< param repoURL >}}/blob/main/tools/render_doc_plots.py):
   add the name to `ALL_RULES` (always) and to `APPLICABLE` (only if
   Layer B passes). Then regenerate:
   ```bash
   python3 tools/render_doc_plots.py
   ```
   The PNGs land under `docs/static/plots/rules/`.

## Step 8 — Build the doc site

```bash
hugo --source docs --minify --destination public
```

Hugo emits to `docs/public/`. Open `docs/public/index.html` in a browser
and click through:

- Home → Rules → your rule's page renders with stencil + convergence plots.
- Home → Families / Grid families → your rule appears under both facets.
- Search the rule name in the address bar — the URL is
  `/rules/<rule_name>/` thanks to the permalink rule in
  [`docs/hugo.toml`]({{< param repoURL >}}/blob/main/docs/hugo.toml).

## Step 9 — Final verification

Run the full test suite once more before opening the PR:

```bash
julia --project=. test/runtests.jl
```

Confirm `test/junit-esd.xml` lists your rule with the outcomes you expect
across all four layers.

## Next steps

- **Walker layer D — discrete conservation.** Conservation acceptance
  (mass / momentum / energy) for finite-volume rules is on the roadmap
  as a parallel layer to B′. The fixture kind will live at
  `fixtures/conservation/` once seeded; track the ESD walker docs and
  `test/walk_esd_tests.jl` for the entry point name when it lands.
- **Cross-binding rule evaluation.** The Julia `eval_coeff` path is one of
  three bindings (Julia, Python, TypeScript) that all consume the same
  AST and must agree byte-for-byte on canonical-form output. The bead
  series **dsc-ve1 / dsc-k86 / dsc-e1g** introduces the cross-binding
  conformance harness. When you add a rule, the byte-equality contract
  applies automatically — no per-binding code is needed because the rule
  is pure declarative AST.
- **Rule × grid matrix.** The rule catalog is expected to grow into a
  rule-by-grid coverage matrix (sibling of the catalog T3 rollup):
  every (rule, grid_family) cell is either green (Layer B passes), amber
  (`applicable: false` with a tracked reason), or red (untested).
  Authoring a new sibling for an additional grid family — e.g.
  `centered_2nd_uniform_<my_family>.json` — is the same nine-step flow
  above, with the per-family selector kind and bindings swapped in.
- **Selector schema design.** If your rule needs a structural concept
  the existing per-family kinds do not express, **read
  [`discretizations/SELECTOR_KINDS.md`]({{< param repoURL >}}/blob/main/discretizations/SELECTOR_KINDS.md)
  before extending the schema.** The decisions log records both the
  resolution and the rejected alternatives for every selector-shape
  question that has come up so far.
