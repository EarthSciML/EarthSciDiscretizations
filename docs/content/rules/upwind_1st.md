---
title: "upwind_1st"
slug: "upwind_1st"
family: "finite_difference"
grid_family: "cartesian"
rule_kind: "scheme"
accuracy: "O(dx)"
applies_to: "grad(u), dim=x"
rule_path: "discretizations/finite_difference/upwind_1st.json"
description: "First-order upwind (backward) finite difference for ∂u/∂x on a uniform Cartesian axis."
tags: ["finite-difference", "upwind", "uniform", "first-order"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/upwind_1st-stencil.png"
       alt="upwind_1st stencil — −1/dx at i−1, +1/dx at i">
  <figcaption>One-sided upwind stencil. The branch shown corresponds to the
  positive-velocity case (information flowing from <code>i−1</code> to <code>i</code>);
  the negative-velocity branch is the mirror image (offsets <code>0</code> and <code>+1</code>
  with sign-flipped coefficients).</figcaption>
</figure>

## Coefficients

| selector kind | axis | offset | coeff |
|---|---|---:|---|
| `cartesian` | `$x` | −1 | `−1 / dx` |
| `cartesian` | `$x` |  0 | `+1 / dx` |

Combined as `+`. Adds an explicit `(dx/2) · u′′(x)` numerical-diffusion term;
this stabilizes advection-dominated problems but is the source of well-known
upwind smearing.

## Convergence

<figure class="figure">
  <img src="/plots/rules/upwind_1st-convergence.png"
       alt="Empirical convergence — slope ≈ −1">
  <figcaption>L∞ error of upwind-1 on <code>u(x) = sin(2πx)</code> with N up to 512.
  The slope cleanly matches the expected −1 reference.</figcaption>
</figure>

Fixture: [`discretizations/finite_difference/upwind_1st/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/upwind_1st/fixtures/convergence)
