---
title: "periodic_bc"
slug: "periodic_bc"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "boundary_condition"
accuracy: "exact (index rewrite)"
applies_to: "index(u, i − offset) on a periodic axis"
rule_path: "discretizations/finite_difference/periodic_bc.json"
description: "Periodic boundary handling as a downstream index-rewrite rule: wraps an out-of-range stencil access index(u, i − offset) to index(u, mod(i − offset + Nx, Nx)) on axes declared periodic by the domain."
tags: ["finite-difference", "boundary-condition", "periodic", "index-rewrite", "P0"]
---

## Rewrite

`periodic_bc` is **not** a discretization scheme — it is one of the
**downstream boundary-condition rewrite rules** that consume the
domain's `boundary_conditions` block (esm-spec §11.5) and rewrite a
lowering's concrete index expressions at the boundary cells. A scheme
rule (e.g. [`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}})
or [`weno5_advection`]({{< ref "/rules/weno5_advection" >}})) emits the
**interior** closed form with bare `index(u, i ± k)` accesses; this rule
fires afterward on every such access whose axis the domain marks
`periodic`, wrapping the index back into range.

The rule matches an `index` access whose index expression is a
metavariable shift `$i − $offset` and replaces it with the modulo-wrapped
form:

```text
applies_to:  index($u, $i − $offset, $rest)
where:       var_has_grid($u, $g)          # $u lives on grid $g
             dim_is_periodic(x, $g)         # axis x of $g is periodic
replacement: index($u, mod(($i − $offset) + Nx, Nx), $rest)
region:      $side
```

`$rest` carries any trailing indices unchanged (so the rule is
axis-local on a multi-dimensional array), and `region: $side` scopes the
rewrite to the boundary side that triggered it. `Nx` is the periodic
axis length supplied by the grid.

After substitution, a stencil neighbor that would read off the low or
high end of the axis instead reads the wrapped-around interior cell:

$$\texttt{index}(u,\; i - \text{off}) \;\longrightarrow\; \texttt{index}\!\left(u,\; (i - \text{off} + N_x) \bmod N_x\right).$$

The `+ Nx` before the `mod` keeps the argument non-negative so the
wrap is correct for left-edge accesses (`i − off < 0`) under
truncating-modulo semantics.

## Guards

Two `where` guards gate the rewrite so it fires only where periodicity
actually applies:

| Guard | Binds / checks |
|---|---|
| `var_has_grid` | binds `$g` to the grid the accessed variable `$u` lives on |
| `dim_is_periodic` | fires only when axis `x` of grid `$g` is declared periodic in the domain's `boundary_conditions` |

Because the guards read the `(grid, axis-periodicity)` pair from the
domain, the same closed scheme lowering composes with a periodic axis
here and with a Dirichlet/Neumann fill row elsewhere — the scheme rule
itself stays BC-agnostic. This is the mechanism referenced by the BC
handoff tables on the scheme pages
([`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}),
[`weno5_advection`]({{< ref "/rules/weno5_advection" >}})).

## Convergence

This is an index-rewrite rule, not a manufactured-solution scheme, so it
has no standalone empirical order of accuracy — periodic wrapping is
*exact* (it relabels indices; it introduces no truncation error). Its
acceptance signature is therefore a **rewrite fixture** rather than an
MMS convergence sweep: the Layer-B convergence fixture under
[`discretizations/finite_difference/periodic_bc/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/periodic_bc/fixtures/convergence)
carries `applicable: false` (with a tracked `skip_reason`), and the
authoritative check is the AST-rewrite golden under
[`fixtures/rewrite/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/periodic_bc/fixtures/rewrite),
which pins `index(T, i − 1, j)` on an `x`-periodic grid rewriting to the
`mod`-wrapped form. The correctness of the *scheme* that this rule
supports is verified on a periodic domain by that scheme's own
convergence fixture (e.g. the `centered_2nd_uniform` sweep on
`sin(2πx)`).
