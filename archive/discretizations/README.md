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
  (also holds the boundary-condition rules: `dirichlet_bc.json`,
  `neumann_bc.json`, `periodic_bc.json`, `robin_bc.json`)
- `finite_volume/` — finite-volume reconstructions and flux forms
- `ic/` — initial-condition rules
- `spectral/` — spectral / pseudo-spectral methods. No catalog rules yet: the
  declarative feasibility of the family is characterized in
  `spectral/SPECTRAL_FEASIBILITY.md` (Fourier/Chebyshev differentiation are
  AST-expressible and proven to machine precision but blocked from being a
  general rule by an ESS reduction-range-sizing gap; spherical harmonics are
  fundamentally infeasible). See that verdict before adding spectral rules.
- `integral/` — quadrature operators (domain integrals / PIDE reduction terms,
  integral operators). Ships `midpoint_1d.json` (full-domain 1-D midpoint
  quadrature as a `reduce` sum_product arrayop). The family's declarative
  feasibility is characterized in `integral/INTEGRAL_FEASIBILITY.md`: the
  reduction tracks the grid size via a host-supplied `size_x` const_array
  (`index(size_x, 1)`), bypassing the reduction-range-sizing gap ("GAP A") that
  blocked the spectral family — the same host-supplied-grid-data pattern
  `nn_diffusion` uses. Proven to machine precision (`integral/feasibility_probe/`).
  The bypass is operator-agnostic, so it also unblocks the spectral family. See
  that verdict before adding integral rules.

The per-family split is a starting convention and may evolve as content lands.
See `../docs/REPO_LAYOUT.md` for the governing convention.

## Adding a rule

1. Drop a `.json` file into the appropriate family subdirectory (+ any fixtures under `<name>/fixtures/`)
2. Ensure it validates against the ESS §7 schema
3. Add a rule-specific test under `../test/` if needed (e.g. `test_rule_catalog.jl` structural checks)

**Do NOT:**
- Commit `../docs/rule-catalog.md` — it is **generated** by `../docs/generate_rule_catalog.jl`
  and is gitignored. Running `julia --project=. docs/generate_rule_catalog.jl` regenerates
  it; do not hand-edit or commit it.
- Edit `../test/test_esd_walker.jl` — the walker discovers rules and fixture files
  **dynamically** via `load_rules`. Adding a rule JSON requires no walker edit.

## How a rule reaches the pipeline

A rule is authored as a **closed `replacement` ExpressionNode AST** in
the ESS §4.2 op vocabulary (`arrayop`, `makearray`, `index`,
`broadcast`, arithmetic, `ifelse`, …). The ESS rule engine (`parse_rule`
/ `discretize`) matches the rule's `applies_to` / `pattern` against the
equation tree and substitutes the `replacement` directly, then lifts it
to an `ArrayOp` and evaluates it via the ESS tree-walk evaluator —
the canonical `discretize → ArrayOp → eval` pipeline. There is no ESD
lowering pass: ESD ships the JSON, ESS applies it.

Multi-axis schemes (e.g. `laplacian($u)`) and multi-output schemes
(e.g. `ppm_reconstruction`, which emits both edge fields) express their
structure inside the `replacement` itself — multiple `arrayop` outputs
for a multi-output rule, indexed `arrayop` bodies for multi-axis stencils
— rather than via any external scheme-expansion step. See
[`SELECTOR_KINDS.md`](SELECTOR_KINDS.md) for the per-family selector
vocabulary and the [authoring tutorial](../docs/content/tutorials/add-a-rule/index.md)
for an end-to-end walkthrough.

> **Historical note.** Earlier revisions of ESD carried a
> `src/stencil_lowering.jl` with `lower_stencil_to_scheme` /
> `lower_stencil_to_canonical_replacement` / `lower_stencil_to_replacement`
> helpers that transformed a compact stencil form into a replacement AST.
> That file and those functions have been **retired** — rules are now
> authored directly in replacement form and applied by the ESS rule
> engine, so no in-ESD lowering step exists.
