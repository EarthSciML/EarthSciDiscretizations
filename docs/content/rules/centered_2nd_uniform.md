---
title: "centered_2nd_uniform"
slug: "centered_2nd_uniform"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx²)"
applies_to: "grad(u), dim=x"
rule_path: "discretizations/finite_difference/centered_2nd_uniform.json"
description: "Two-point centered second-order finite difference for ∂u/∂x on a uniform Cartesian axis, expressed as a closed §4.2 arrayop lowering."
tags: ["finite-difference", "centered", "uniform", "P0"]
---

## Lowering

The rule is the canonical exemplar of an ESD linear scheme that lowers a
§4.2 PDE operator to a **closed `arrayop` expression** in the §4.2 op
vocabulary. There are no scheme-specific coefficient blobs and no
off-spec selector / offset fields — `applies_to` matches against `grad`
(itself a §4.2 op), and the lowering is a single `arrayop` whose body
combines `index`, `+`, `-`, `*`, and `/`. ESS bindings (any binding) can
evaluate the lowering by walking the AST.

```text
applies_to:  grad(u, dim=x)
replacement: arrayop(
               output_idx = [$x],
               expr       = (index($u, $x+1) − index($u, $x−1)) / (2·dx),
               args       = [$u]
             )
```

After pattern-variable substitution (`$u → u`, `$x → i`), the lowering at
each interior cell reduces to

$$\left(\frac{\partial u}{\partial x}\right)_i \;\approx\; \frac{u_{i+1} - u_{i-1}}{2\,\Delta x}.$$

## Boundary conditions

Boundary handling is read from the domain's `boundary_conditions` block
(esm-spec §11.5) at lowering time and applied as **downstream rewrite
rules over concrete indices** — there is no `bc:*` op embedded in the
lowered AST. The pattern is:

| Domain BC | Index transformation applied to `$u[$x ± 1]` |
|---|---|
| `periodic` | `mod($x ± 1 + N, N)` (see [`periodic_bc`]({{< ref "/rules/periodic_bc" >}})) |
| `dirichlet` / `constant` | Boundary cells read the prescribed value via `index` into a fill row |
| `zero_gradient` / `neumann` | Mirror the in-range neighbor (`min`/`max` clamp on the index) |

Each transformation is a separate rule that fires at the boundary cells.
The `centered_2nd_uniform` rule itself stays BC-agnostic — its
replacement is the interior closed form. The lowering pipeline (rule
application + BC rewrites) takes the `(grid_family, BC list)` pair from
the domain and emits index expressions that respect the declared BCs.

## Truncation derivation

Symmetric Taylor expansion of the two neighbors about \(x_i\),

$$u_{i\pm 1} \;=\; u_i \;\pm\; \Delta x\,u'_i \;+\; \tfrac{\Delta x^{2}}{2}\,u''_i \;\pm\; \tfrac{\Delta x^{3}}{6}\,u'''_i \;+\; \tfrac{\Delta x^{4}}{24}\,u^{(4)}_i \;\pm\; \cdots,$$

cancels the even-order terms when subtracted, giving

$$\frac{u_{i+1} - u_{i-1}}{2\,\Delta x} \;=\; u'_i \;+\; \tfrac{\Delta x^{2}}{6}\,u'''_i \;+\; O(\Delta x^{4}).$$

The scheme is therefore second-order accurate, with a purely dispersive
leading error \(\tfrac{\Delta x^{2}}{6}\,u'''(x)\) and no numerical
diffusion — contrast with [`upwind_1st`]({{< ref "/rules/upwind_1st" >}}),
whose one-sided stencil introduces an explicit diffusive
\((\Delta x / 2)\,u''(x)\) term.

## Convergence

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform-convergence.png"
       alt="Empirical convergence — slope ≈ −2 on log-log">
  <figcaption>L∞ error of the centered stencil applied to <code>u(x) = sin(2πx)</code>
  on a periodic <code>[0, 1]</code> domain, sampled at cell centers. Empirical slope
  matches the expected −2 reference line.</figcaption>
</figure>

The fixture under
[`discretizations/finite_difference/centered_2nd_uniform/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform/fixtures/convergence)
sets `expected_min_order = 1.9` to tolerate minor pre-asymptotic drift on
the 16 → 32 → 64 → 128 sequence. The sweep is currently marked
`applicable: false` pending an ESS upgrade that lets `mms_convergence`
walk a `replacement`-form rule through the existing arrayop / broadcast
evaluator (no new scheme-specific kernels). Once ESS lands the AST
dispatch path, the convergence fixture re-enables without modification.
