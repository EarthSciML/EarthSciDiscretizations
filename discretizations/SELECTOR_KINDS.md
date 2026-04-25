# Selector Kinds — Rule Schema Index

This document is the running index of stencil-selector design decisions made
while porting discretization rules across grid families. Each row records
**one** resolved question, the chosen convention, and the bead that pinned it.
New rule authors should read this before introducing a new grid family.

The selector schema lives inside the ESS §7 rule format:

```json
{
  "stencil": [
    {
      "selector": { "kind": "<family>", "axis": "<axis>", "offset": <int> },
      "coeff":    { "...AST..." }
    }
  ]
}
```

`kind` binds the offset to a specific grid-family accessor. `axis` names the
direction along which `offset` is interpreted. `coeff` is an
`ExpressionNode` AST evaluated against per-grid bindings (see ESS §9.2).

## Decisions

| # | Question | Decision | Pinned by |
|---|---|---|---|
| 1 | Selector `kind` naming: per-family (`vertical`, `latlon`, `cubed_sphere`) or structural (`structured_1d`, `structured_2d_metric`)? | **Per-family.** Each kind names exactly one grid-family accessor; structural categories may emerge once 3+ families share a stencil but are premature today. | dsc-mzu |
| 2 | Spacing symbol for uniform 1D structured grids (cartesian, vertical-uniform). | **`h`** is the canonical generic uniform-spacing symbol. Per-axis aliases (`dx`, `dz`, …) are also bound by the MMS harness so existing cartesian rules keep working. New rules SHOULD prefer `h` to make grid-family-agnostic AST reuse possible. | dsc-mzu |
| 3 | Metric for non-uniform vertical grids (`eta`, `theta`, `z` stretched). | **Defer.** The convention will be a per-cell symbolic reference `dz_k` resolved against the grid accessor at evaluation time. Implementing it requires extending the ESS evaluator + grid binding contract; tracked as a follow-up bead. Until then, only uniform vertical (`sigma_uniform`, `z_uniform`) is supported by Layer B. | dsc-mzu (deferred to follow-up) |
| 4 | Boundary handling for non-periodic vertical (top/bottom). | **Compose with a separate BC rule** (parallel to `periodic_bc.json`) rather than inlining BC handling in the stencil. Keeps the centered stencil rule clean and lets the same BC rule cover all 1D structured families. Concrete `vertical_top_bottom_bc.json` is a follow-up. | dsc-mzu (deferred to follow-up) |
| 5 | One rule with `grid_family: ["cartesian", "vertical"]` vs sibling files. | **Sibling files** when the selector `kind` differs (which is always, given decision #1). A `grid_family` array is reserved for cases where the AST and selectors are byte-identical — currently no such case exists. | dsc-mzu |

## When to add a new selector kind

Adding a new family (e.g. `latlon`, `cubed_sphere`):

1. Pick a per-family kind name in lowercase snake_case.
2. Decide the canonical axis names (`(:lon, :lat)`, `(:xi, :eta)`, …).
3. Pick the spacing symbol convention. Reuse `h` for axis-agnostic uniform
   spacing; introduce per-axis symbols (`dlon`, `dlat`) only when metrics
   require it.
4. Add a row above with the rationale and the bead that pinned it.
5. Author the first rule in `discretizations/<family>/` with a paired
   convergence fixture under `<rule>/fixtures/convergence/`.

## See also

- `../docs/rule-catalog.md` — manifest of rules to port across families.
- `README.md` (this directory) — rule-file authoring policy.
- ESS `esm-spec.md` §7 — rule schema; §9.2 — `call`-op decision tree.
