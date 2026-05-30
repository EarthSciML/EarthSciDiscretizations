# Agent instructions — EarthSciDiscretizations

This file documents project-wide expectations for any agent authoring code
or discretization rules in this repository. It is the committed counterpart
to the per-agent scratch `CLAUDE.md` (which is gitignored and does not ship
to future sessions).

## Authoring discretization rules

Rules live in `discretizations/<family>/*.json`. The schema is ESS §7
(see `EarthSciSerialization/esm-spec.md`). **Express all math directly
in the ExpressionNode AST** — available ops include `+ - * / ^`,
`max min ifelse abs sign`, all trig/exp/log, comparisons, and logical
ops. See ESS `esm-spec.md` §9.2 decision tree for when (rarely) a
`call` op is justified.

Do NOT:

- Write Julia/Python/etc helper functions that implement math a rule
  should own (limiters, reconstructions, ratio functions).
- Use `{"op": "call", "fn": "…"}` for anything expressible in existing
  ops.
- Rely on per-binding runtime names that hide the math.

When in doubt: if the formula fits on paper, it fits in the AST.

See `discretizations/README.md` for the catalog landing and the
[generated rule catalog](docs/rule-catalog.md) for the full manifest.

### What rule authors must NOT commit

- `docs/rule-catalog.md` — **generated** by `docs/generate_rule_catalog.jl` at
  doc-build time; gitignored. Do not commit it; do not run the generator and
  stage its output. The CI doc build regenerates it automatically.
- Edits to `test/test_esd_walker.jl` — the walker discovers rules and fixture
  files **dynamically**. Adding a rule JSON (with or without fixtures) requires
  zero edits to the walker test.

When authoring a new rule, the ONLY files you commit are:
1. `discretizations/<family>/<name>.json` — the rule itself
2. `discretizations/<family>/<name>/fixtures/**` — any fixtures (canonical,
   convergence, rewrite, monotonicity, conservation)
3. Rule-specific test items in `test/test_rule_catalog.jl` (structural checks)

## The single-pathway rule (all bindings)

**ESD never carries a rule evaluator. In every binding, the rule
evaluator/runner is a THIN PASSTHROUGH to the corresponding ESS
binding's evaluator.** There is exactly one canonical pipeline:

> rule application in ESS → ArrayOp → eval

No binding may carry a shadow evaluator, a reimplementation of an ESS
op, or a binding-local fast path that bypasses ESS. If you find one,
file a bead to retire it; do not extend it.

Per-binding entry points (each is a passthrough to ESS):

| Binding    | ESD entry point                                           | Delegates to                                                      |
|------------|-----------------------------------------------------------|-------------------------------------------------------------------|
| Julia      | `EarthSciDiscretizations.eval_coeff` (`src/rule_eval.jl`) | `EarthSciSerialization.evaluate`                                  |
| Python     | `earthsci_discretizations.rules` (`python/src/earthsci_discretizations/rules/`) | the `earthsci_serialization` Python evaluator        |
| Rust       | `rule_eval` (`rust/src/rule_eval.rs`)                     | the `earthsci_serialization` Rust evaluator                       |
| TypeScript | `rules/` (`typescript/src/rules/`)                        | the `earthsci_serialization` TypeScript evaluator                 |

A binding-level file with non-trivial math, branch logic, or op
dispatch is a bug, not a feature. Passthrough means: marshal inputs,
call ESS, return outputs. Nothing else.

### Conformance: golden regeneration drives the canonical pipeline

The `regenerate_golden.*` scripts under `tests/conformance/` (e.g.
`tests/conformance/grids/cartesian/regenerate_golden.jl`,
`tests/conformance/grids/latlon/regenerate_golden.py`,
`tests/conformance/rules/centered_2nd_uniform_latlon/regenerate_golden.jl`)
**MUST drive the canonical pipeline** — rule application in ESS →
ArrayOp → eval through the ESD passthrough. They must NOT call a
per-binding shadow evaluator. If a regeneration script computes
expected values via a binding-local code path, the resulting golden
will mask divergence between bindings rather than expose it.

### Cross-reference

The single-pathway rule is the project-level expression of the
workspace-root `CLAUDE.md` (Gas Town mayor) directive that ESD is a
discretization catalog over ESS, not a parallel runtime.

## Dependency resolution

`EarthSciSerialization` is not yet in the Julia General registry. ESD's
`Project.toml` declares a Pkg `[sources]` entry pointing at the GitHub
repo (`packages/EarthSciSerialization.jl` subdir on `main`), so a plain
`Pkg.instantiate()` resolves ESS directly from GitHub on Julia ≥ 1.11.
The project's `julia` compat is therefore `"1.11"`.

For local development against an in-progress ESS checkout (Gas Town
workspace), `scripts/setup_polecat_env.sh` Pkg.develops a local path —
the explicit develop entry overrides the `[sources]` URL until you
remove it. The script is no longer required for CI.
