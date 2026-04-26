---
title: "weno5_advection"
slug: "weno5_advection"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
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

Numeric coverage lives in the canonical Julia test fixture at
[<code>tests/fixtures/weno5_advection/</code>]({{< param repoURL >}}/blob/main/tests/fixtures/weno5_advection),
exercised by
[<code>test/test_weno5_advection_rule.jl</code>]({{< param repoURL >}}/blob/main/test/test_weno5_advection_rule.jl).
The MMS fixture uses `f(x) = sin(2πx + 1)` on `[0, 1]` with periodic boundary
conditions; the phase shift keeps the critical points of `f` away from every
dyadic cell face at `n ∈ {32, 64, 128, 256}`, sidestepping the well-known
WENO5-JS accuracy dip from `ω_k → d_k` recovery stalling at
`f'(x_{i+1/2}) = 0` (Henrick, Aslam & Powers, JCP 2005).

| `n`  | `dx`         | L∞ error     | observed order |
| ---: | :----------- | :----------- | :------------- |
|  32  | 0.03125000   | 3.4137e-05   | —              |
|  64  | 0.01562500   | 1.0656e-06   | 5.002          |
| 128  | 0.00781250   | 3.3254e-08   | 5.002          |
| 256  | 0.00390625   | 1.0377e-09   | 5.002          |

Theoretical asymptotic order: **5.0** (Jiang & Shu 1996, smooth regions).
Acceptance threshold: **min(observed order) ≥ 4.7** — leaves headroom for the
small accuracy hit from the `ε = 1e-6` regularisation of the nonlinear
weights. A companion shock-capturing fixture advects a unit square wave one
full period at CFL 0.4 with SSP-RK3; max overshoot/undershoot is ~3.7e-4,
well under the 0.05 tolerance.

<div class="callout callout-pending">
<strong>Pending ESS walker harness extension.</strong>
The Layer-B walker harness needs nonlinear-reconstruction support — the
ratio-form nonlinear weights `ω_k = d_k / (ε + β_k)²` are not yet expressible
in the §7 stencil schema. The walker-side fixture under
[<code>discretizations/finite_volume/weno5_advection/fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/weno5_advection/fixtures/convergence)
records a structured SKIP and will populate once that extension lands; until
then the canonical numeric coverage above is the source of truth.
</div>
