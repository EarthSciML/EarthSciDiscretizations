---
title: "centered_2nd_uniform_vertical"
slug: "centered_2nd_uniform_vertical"
families: "finite_difference"
grid_families: "vertical"
rule_kinds: "scheme"
accuracy: "O(h²)"
applies_to: "grad(u), dim=k"
rule_path: "discretizations/finite_difference/centered_2nd_uniform_vertical.json"
description: "Two-point centered second-order finite difference for ∂u/∂k on a uniformly-spaced vertical axis."
tags: ["finite-difference", "centered", "vertical", "uniform"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform_vertical-stencil.png"
       alt="centered_2nd_uniform_vertical stencil">
  <figcaption>Two-point centered stencil along the vertical axis.</figcaption>
</figure>

## Coefficients

| selector kind | axis | stagger | offset | coeff |
|---|---|---|---:|---|
| `vertical` | `$k` | `face_bottom` | 0 | `−1 / h` |
| `vertical` | `$k` | `face_top`    | 0 | `+1 / h` |

Combined as `+`. The `vertical` selector kind addresses face-staggered
samples of `u`: at cell `i`, `face_bottom` with offset `0` reads the face
at the bottom of the cell, and `face_top` with offset `0` reads the face at
the top of the cell. The vertical accessor exposes the level spacing `h`
directly so the rule does not need to know whether the column is in `z`,
`p`, or `σ` coordinates.

## Discrete operator

For face samples $u_{i+1/2}$ on a uniform vertical axis with cell spacing
$h$, the rule produces a cell-centered derivative at cell $i$ (located at
$z_i$, the midpoint of the cell's two faces $z_{i-1/2}$ and $z_{i+1/2}$):

$$\bigl(\partial_k u\bigr)_i \;\approx\; \frac{u_{i+1/2} - u_{i-1/2}}{h}.$$

Taylor-expanding $u_{i\pm 1/2}$ about cell center $z_i$,

$$\frac{u_{i+1/2} - u_{i-1/2}}{h} \;=\; u'(z_i) \;+\; \frac{h^{2}}{24}\,u'''(z_i) \;+\; \mathcal{O}(h^{4}),$$

so the leading truncation error is $\tfrac{h^{2}}{24}\,u'''(z_i)$, giving
the advertised $\mathcal{O}(h^{2})$ accuracy. The expression is
anti-symmetric in the two adjacent faces, so the operator is
non-dissipative — error appears as dispersion only, with no built-in
numerical diffusion.

The vertical accessor binds `h` to whichever coordinate the column actually
uses ($\Delta z$, $\Delta p$, or $\Delta\sigma$); the discrete operator
above is unchanged. Non-uniform spacing is out of scope for this rule —
the uniform axis is enforced by the `vertical` selector kind.

## Convergence

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform_vertical-convergence.png"
       alt="Empirical convergence on vertical axis">
  <figcaption>L∞ error on <code>u(z) = sin(2πz)</code>, column <code>[0, 1]</code>,
  face-staggered samples → cell-centered derivative. Slope matches the −2 reference.</figcaption>
</figure>

Fixture: [`discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence)
