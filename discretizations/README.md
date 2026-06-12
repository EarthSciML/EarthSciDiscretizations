# Discretization Rule Files

This directory holds authoritative discretization rule JSON files. Each file
describes how a continuous PDE operator maps onto a discrete stencil, grid
staggering, and boundary-condition handling.

Rule files are validated against the EarthSciSerialization discretization
schema (§7). Validation and application are performed by the live ESS rule
engine: `parse_rule` consumes replacement-form rules, coefficient ASTs
evaluate via `eval_coeff` (a passthrough to ESS
`parse_expression`/`evaluate_expr`), and `discretize → ArrayOp → eval`
drives the canonical pipeline exercised by `../test/walk_esd_tests.jl`.

**AST-first authoring policy.** Express all math directly in the
ExpressionNode AST (ESS §7; see `../AGENTS.md` "Authoring discretization
rules" for the full policy). Coefficient and reconstruction terms
evaluate via `EarthSciDiscretizations.eval_coeff`, a thin passthrough to
the ESS tree-walk evaluator — ESD does not carry a shadow evaluator. If
the formula fits on paper, it fits in the AST; reach for `{op: "call"}`
only when the ESS `esm-spec.md` §9.2 decision tree says so.

## Layout

- `finite_difference/` — finite-difference stencils on structured grids
- `finite_volume/` — finite-volume reconstructions and flux forms
- `spectral/` — spectral / pseudo-spectral methods

The per-family split is a starting convention and may evolve as content lands.
See `../docs/REPO_LAYOUT.md` for the governing convention.

## Adding a rule

1. Drop a `.json` file into the appropriate family subdirectory (+ any fixtures under `<name>/fixtures/`)
2. Ensure it validates against the ESS §7 schema
3. Add a rule-specific test under `../test/` if needed (e.g. `test_rule_catalog.jl` structural checks)

**Do NOT:**
- Commit `../docs/rule-catalog.md` — it is **generated** by `../docs/generate_rule_catalog.jl`
  at doc-build time and is gitignored. Running `julia --project=docs docs/make.jl` regenerates
  it; do not hand-edit or commit it.
- Edit `../test/test_esd_walker.jl` — the walker discovers rules and fixture files
  **dynamically** via `load_rules`. Adding a rule JSON requires no walker edit.

## Stencil-form vs replacement-form rules

Rules are authored in one of two ESS §7 forms:

- **Replacement form** — a single `replacement` ExpressionNode AST that
  ESS `rule_engine` (`parse_rule`) consumes directly. Used by rules whose
  full math fits compactly inline (e.g. `centered_2nd_uniform.json`).
- **Stencil form** — an array of `{selector, coeff}` entries plus a
  `combine` op. Compact and grid-family-aware (see
  [`SELECTOR_KINDS.md`](SELECTOR_KINDS.md)). Used by upwind / flux /
  reconstruction rules whose shape varies by stagger (e.g.
  `upwind_1st.json`, `lax_friedrichs_flux.json`).

A stencil-form rule reaches the ESS pipeline by one of three lowerings:

- [`EarthSciDiscretizations.lower_stencil_to_scheme`](../src/stencil_lowering.jl)
  — emits the ESM document parts (`discretizations.<name>` scheme +
  `use:` rule, RFC §7.2.1) that ESS `discretize` consumes natively via
  scheme expansion → ArrayOp lift → tree-walk eval. This is the
  canonical-pipeline path; cartesian-only today, mirroring the ESS
  scheme-expansion foundation (esm-j1u).
- [`EarthSciDiscretizations.lower_stencil_to_canonical_replacement`](../src/stencil_lowering.jl)
  (dsc-vst2) — a replacement expression whose axis references use ESS's
  canonical `$target` component names (`i`, `j`, …; RFC §7.1,
  position-based). The form multi-axis rules with dim-less patterns
  (e.g. `laplacian($u)`) need to ride ESS `discretize` as plain
  `pattern` + `replacement` document rules; collocated (cell_center)
  stencils only.
- [`EarthSciDiscretizations.lower_stencil_to_replacement`](../src/stencil_lowering.jl)
  (dsc-y0jj) — inserts a scalar `replacement` AST into the rule dict
  (legacy path; supports more selector families, but downstream
  consumers must materialize indices themselves). Currently supports the
  `cartesian`, `arakawa`, `latlon`, `cubed_sphere`, and `vertical`
  selector kinds; other kinds (`indirect`, `reduction`, …) raise
  `ArgumentError` until their lowering rows are added.

Both lowerers are pure AST → AST transforms — no per-rule-shape
dispatch, no numerical evaluation; each new selector kind composes as a
separate dispatch branch.
