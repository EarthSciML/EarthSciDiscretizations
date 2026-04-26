---
title: "centered_2nd_uniform_vertical"
slug: "centered_2nd_uniform_vertical"
family: "finite_difference"
grid_family: "vertical"
rule_kind: "scheme"
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

| selector kind | axis | offset | coeff |
|---|---|---:|---|
| `vertical` | `$k` | −1 | `−1 / (2 h)` |
| `vertical` | `$k` | +1 | `+1 / (2 h)` |

Combined as `+`. Same kernel as
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}),
relabeled for the vertical axis. The vertical accessor exposes the level
spacing `h` directly so the rule does not need to know whether the column
is in `z`, `p`, or `σ` coordinates.

## Discrete operator

For a column sample $u_k$ on a uniform vertical axis with spacing $h$, the
rule produces

$$\bigl(\partial_k u\bigr)_k \;\approx\; \frac{u_{k+1} - u_{k-1}}{2\,h}.$$

Taylor-expanding $u_{k\pm 1}$ about level $k$,

$$\frac{u_{k+1} - u_{k-1}}{2\,h} \;=\; u'_k \;+\; \frac{h^{2}}{6}\,u'''_k \;+\; \mathcal{O}(h^{4}),$$

so the leading truncation error is $\tfrac{h^{2}}{6}\,u'''(k)$, giving the
advertised $\mathcal{O}(h^{2})$ accuracy. The expression is anti-symmetric
in $k\!\pm\!1$, so the operator is non-dissipative — error appears as
dispersion only, with no built-in numerical diffusion.

The vertical accessor binds `h` to whichever coordinate the column actually
uses ($\Delta z$, $\Delta p$, or $\Delta\sigma$); the discrete operator
above is unchanged. Non-uniform spacing is out of scope for this rule —
the uniform axis is enforced by the `vertical` selector kind.

## Convergence

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform_vertical-convergence.png"
       alt="Empirical convergence on vertical axis">
  <figcaption>L∞ error on <code>u(z) = cos(2πz)</code>, periodic <code>[0, 1]</code>,
  cell-centered samples. Slope matches the −2 reference.</figcaption>
</figure>

Fixture: [`discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence)
