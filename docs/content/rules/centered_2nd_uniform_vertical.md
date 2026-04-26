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

Same kernel as [`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}),
relabeled for the vertical axis. The vertical accessor exposes the level
spacing `h` directly so the rule does not need to know whether the column
is in `z`, `p`, or `σ` coordinates.

## Convergence

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform_vertical-convergence.png"
       alt="Empirical convergence on vertical axis">
  <figcaption>L∞ error on <code>u(z) = cos(2πz)</code>, periodic <code>[0, 1]</code>,
  cell-centered samples. Slope matches the −2 reference.</figcaption>
</figure>

Fixture: [`discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_vertical/fixtures/convergence)
