---
title: "Discretization rules"
description: "Closed §4.2 AST lowerings for the seeded ESD rule files, with boundary handling delegated to the domain."
---

The seeded rule files. Each rule **lowers a §4.2 PDE operator** —
`grad`, `div`, `laplacian`, `D`, plus pointwise math — to a **closed
expression in the §4.2 op vocabulary**: `arrayop`, `broadcast`,
`ifelse`, `+`, `-`, `*`, `/`, `^`, `sqrt`, and friends. There are no
scheme-specific kernels in any binding language, and `applies_to.op`
only matches §4.2 ops (canonical names like `grad` and `div`; never
off-spec selectors such as `advect`, `reconstruct`, `flux`, or
`limit`). Any ESS binding can evaluate a lowering by walking its AST.

## Boundary handling

Boundary conditions live on the domain, not in the rule. The domain's
`boundary_conditions` block (esm-spec §11.5) declares the BC type per
axis — `periodic`, `dirichlet`, `neumann`, `zero_gradient`, or
`robin` — and downstream BC rewrite rules (e.g.
[`periodic_bc`]({{< ref "/rules/periodic_bc" >}})) consume that list
and rewrite the lowered AST into concrete index expressions at the
boundary cells. The lowering itself contains no `bc:*` nodes; it is
the interior closed form, and the `(grid_family, BC list)` pair from
the domain drives the boundary rewrites.

## Migration status

The catalog is mid-migration to this closed-AST authoring pattern.
Roughly 14 of the 16 catalog rules still use the legacy
stencil/coefficient form (an explicit `stencil` block plus per-offset
`coefficients`) and are scheduled for rewrite. Their pages describe
the legacy form until the rule itself migrates.

- [`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}})
  is the **canonical linear exemplar** — a single `arrayop` whose body
  combines `index`, `+`, `-`, `*`, and `/` to express
  $(u_{i+1} - u_{i-1}) / (2\,\Delta x)$.
- [`weno5_advection`]({{< ref "/rules/weno5_advection" >}}) will be
  the **canonical nonlinear exemplar** once its rewrite lands (blocked
  on ESS work tracked under `esm-4gw`).

## Per-page rendering

After migration, each rule page renders the closed lowering AST
diagrammatically — pattern-variable bindings on `applies_to`, then the
op tree of the `replacement` — rather than a coefficient diagram.
Pages for rules that have not yet migrated continue to describe the
legacy stencil/coefficient form; they update when their rule does.
Convergence plots (when present) report empirical order of accuracy on
a manufactured solution; rules whose Layer-B fixtures depend on
in-flight ESS harness extensions show a *pending* placeholder until
those fixtures land.
