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
  should own (limiters, reconstructions, ratio functions). The single
  supported evaluator entry point in ESD is
  `EarthSciDiscretizations.eval_coeff`, a thin passthrough to the
  EarthSciSerialization tree-walk evaluator
  (`EarthSciSerialization.evaluate`). ESD does not carry a shadow
  evaluator.
- Use `{"op": "call", "fn": "…"}` for anything expressible in existing
  ops.
- Rely on per-binding runtime names that hide the math.

When in doubt: if the formula fits on paper, it fits in the AST.

See `discretizations/README.md` for the catalog landing and
`docs/rule-catalog.md` for the full rule manifest. The evaluator
passthrough source is `src/rule_eval.jl`.

## Dependency resolution

`EarthSciSerialization` is not yet in the Julia General registry. Both
local development and CI resolve it via `scripts/setup_polecat_env.sh`,
which prefers a Gas Town workspace checkout and falls back to
`Pkg.add(url=...)` from GitHub. Run the script once before
`julia --project=. -e 'using EarthSciDiscretizations'` in a fresh
environment.
