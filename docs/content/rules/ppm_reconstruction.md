---
title: "ppm_reconstruction"
slug: "ppm_reconstruction"
family: "finite_volume"
grid_family: "cartesian"
rule_kind: "scheme"
accuracy: "O(dx³)"
applies_to: "reconstruct(q), dim=x"
rule_path: "discretizations/finite_volume/ppm_reconstruction.json"
description: "Piecewise-parabolic method (Colella & Woodward 1984) — 4th-order edge values + a parabola per cell, no limiter."
tags: ["finite-volume", "ppm", "reconstruction", "colella-woodward"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/ppm_reconstruction-stencil.png"
       alt="PPM 4-point edge-value stencil">
  <figcaption>Colella & Woodward (1984) eq. (1.6): the 4th-order edge value
  <code>q<sub>i+1/2</sub></code> is a centered combination of cell averages at
  offsets <code>−1, 0, +1, +2</code>. Two more cells (<code>−2</code> and
  <code>+2</code>) appear in the rule's broader stencil so the same
  edge-value formula applies at <code>q<sub>i−1/2</sub></code>.</figcaption>
</figure>

## Coefficients

Edge value at `q_{i+1/2}`:

| selector | offset | coeff |
|---|---:|---|
| `cartesian` | −1 | `−1/12` |
| `cartesian` |  0 | `+7/12` |
| `cartesian` | +1 | `+7/12` |
| `cartesian` | +2 | `−1/12` |

Per-cell parabola (CW84 eqs. 1.5, 1.7, 1.10):

- `a_L = q_{i−1/2}`  (right limit of cell `i−1`)
- `a_R = q_{i+1/2}`  (left limit of cell `i+1`)
- `da  = a_R − a_L`
- `a₆  = 6·(q_i − ½(a_L + a_R))`
- `a(ξ) = a_L + ξ·(da + a₆·(1 − ξ))`,  `ξ ∈ [0, 1]`

Limiting (CW84 eqs. 1.10) is **not** applied at this rule level — see the
flux-limiter rules ([`flux_limiter_minmod`]({{< ref "/rules/flux_limiter_minmod" >}}),
[`flux_limiter_superbee`]({{< ref "/rules/flux_limiter_superbee" >}})).

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B harness needs sub-stencil targeting (so the edge-value 4-point
stencil and the parabola pass can be exercised independently) plus a parabola
evaluation pass to reconstruct cell-interior values for the manufactured
solution. The fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_reconstruction/fixtures/convergence)
will populate once those extensions land.
</div>
