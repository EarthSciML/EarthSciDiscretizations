---
title: "weno5_advection"
slug: "weno5_advection"
family: "finite_volume"
grid_family: "cartesian"
rule_kind: "scheme"
accuracy: "O(dx⁵)"
applies_to: "advect(q, u), dim=x"
rule_path: "discretizations/finite_volume/weno5_advection.json"
description: "Jiang–Shu 5th-order WENO reconstruction — three sub-stencils, smoothness indicators, nonlinear weights."
tags: ["finite-volume", "weno", "jiang-shu", "advection", "high-order"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/weno5_advection-stencil.png"
       alt="WENO5 three sub-stencils, left-biased branch">
  <figcaption>Three candidate 3-cell sub-stencils for the left-biased
  reconstruction of the edge value <code>q<sub>i+1/2</sub><sup>L</sup></code>.
  Linear weights <code>(d₀, d₁, d₂) = (1/10, 6/10, 3/10)</code> recover formal
  5th-order accuracy in smooth regions; nonlinear weights ω<sub>k</sub> shift
  away from sub-stencils that straddle a discontinuity.</figcaption>
</figure>

## Coefficients (left-biased)

| sub-stencil | offsets | coefficients |
|---|---|---|
| `p₀` | (−2, −1, 0) | (1/3, −7/6, 11/6) |
| `p₁` | (−1, 0, +1) | (−1/6, 5/6, 1/3) |
| `p₂` | (0, +1, +2) | (1/3, 5/6, −1/6) |

Optimal linear weights `(d₀, d₁, d₂) = (1/10, 6/10, 3/10)`. Smoothness
indicators (Jiang & Shu 1996 eq. 2.17):

$$
\beta_0 = \tfrac{13}{12}(q_{i-2}-2q_{i-1}+q_i)^2 + \tfrac14(q_{i-2}-4q_{i-1}+3q_i)^2,
$$

$$
\beta_1 = \tfrac{13}{12}(q_{i-1}-2q_i+q_{i+1})^2 + \tfrac14(q_{i-1}-q_{i+1})^2,
$$

$$
\beta_2 = \tfrac{13}{12}(q_i-2q_{i+1}+q_{i+2})^2 + \tfrac14(3q_i-4q_{i+1}+q_{i+2})^2.
$$

Nonlinear weights `α_k = d_k / (ε + β_k)²`, normalized to
`ω_k = α_k / Σ α_j`. The right-biased branch (used when `u < 0`) is the
mirror of the left-biased branch under index reflection.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B harness needs nonlinear-reconstruction support — the
ratio-form nonlinear weights are not yet expressible in the §7 stencil
schema. The fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/weno5_advection/fixtures/convergence)
will populate once that extension lands.
</div>
