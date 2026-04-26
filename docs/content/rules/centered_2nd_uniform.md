---
title: "centered_2nd_uniform"
slug: "centered_2nd_uniform"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx²)"
applies_to: "grad(u), dim=x"
rule_path: "discretizations/finite_difference/centered_2nd_uniform.json"
description: "Two-point centered second-order finite difference for ∂u/∂x on a uniform Cartesian axis."
tags: ["finite-difference", "centered", "uniform", "P0"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform-stencil.png"
       alt="centered_2nd_uniform stencil — −1/(2dx) at i−1, +1/(2dx) at i+1">
  <figcaption>Two-point centered stencil — neighbors at offset ±1 with
  symmetric coefficients.</figcaption>
</figure>

## Coefficients

| selector kind | axis | offset | coeff |
|---|---|---:|---|
| `cartesian` | `$x` | −1 | `−1 / (2 dx)` |
| `cartesian` | `$x` | +1 | `+1 / (2 dx)` |

Combined as `+`. The result is the cell-centered approximation of the first
spatial derivative with a leading error term `(dx² / 6) · u′′′(x)`.

## Discrete operator

Applied to a cell-centered field \(u_i = u(x_i)\) on a uniform axis with
spacing \(\Delta x\), the rule produces the discrete first-derivative
operator

$$\left(\frac{\partial u}{\partial x}\right)_i \;\approx\; \frac{u_{i+1} - u_{i-1}}{2\,\Delta x}.$$

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
the 16 → 32 → 64 → 128 sequence.
