---
title: "weno5_advection_2d"
slug: "weno5_advection_2d"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(h⁵)"
applies_to: "advect(q, U), dims=(x, y)"
rule_path: "discretizations/finite_volume/weno5_advection_2d.json"
description: "Dimension-by-dimension Jiang–Shu WENO5 reconstruction on a 2D structured uniform Cartesian grid."
tags: ["finite-volume", "weno", "jiang-shu", "advection", "high-order", "two-dimensional"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/weno5_advection_2d-stencil.png"
       alt="2D WENO5 cross stencil: five cells along x and four extra along y, intersecting at (i,j)">
  <figcaption>Dimension-by-dimension splitting (Shu 1998 §2.2; LeVeque 2002
  §20): the 1D Jiang–Shu (1996) WENO5 reconstruction is applied independently
  along the x-axis (blue) and the y-axis (red) through the same cell
  <code>(i,j)</code>. Each axis carries the three candidate sub-stencils of
  the 1D rule. The union spans nine cells per axis (five own-axis +
  four cross-axis), 17 of the 5×5 neighborhood total.</figcaption>
</figure>

## Coefficients

Per-axis blocks under <code>axes.x</code> / <code>axes.y</code> reuse the 1D
[<code>weno5_advection</code>]({{< relref "weno5_advection.md" >}})
sub-stencils, linear weights `(d₀, d₁, d₂) = (1/10, 6/10, 3/10)`, and
smoothness indicators (Jiang & Shu 1996 eq. 2.17), with the selector axis
swapped between `$x` and `$y`. Nonlinear weights `α_k = d_k / (ε + β_k)²`,
normalized to `ω_k = α_k / Σ α_j`. The right-biased branch (used when the
local face velocity is negative) is the index-reflected mirror of the
left-biased branch.

The flux selection per axis is upwinded against the local face velocity:

$$
F^x_{i+1/2,j} = u_{i+1/2,j}\, q_{i+1/2,j}^{\mathrm{upwind}},\qquad
F^y_{i,j+1/2} = v_{i,j+1/2}\, q_{i,j+1/2}^{\mathrm{upwind}}.
$$

The advective tendency sums the two axis flux divergences:

$$
\mathrm{advect}(q, U) =
\frac{F^x_{i+1/2,j} - F^x_{i-1/2,j}}{\Delta x}
\,+\,
\frac{F^y_{i,j+1/2} - F^y_{i,j-1/2}}{\Delta y}.
$$

## Convergence

<figure class="figure">
  <img src="/plots/rules/weno5_advection_2d-convergence.png"
       alt="Empirical L∞ convergence of dimension-by-dimension WENO5 on the 2D phase-shifted sine">
  <figcaption>Cell-averaged inputs of
  <code>u(x,y) = sin(2π x + 1)·sin(2π y + 1/2)</code> on
  <code>[0,1]²</code> periodic, reconstructed face values vs the analytic
  cross-section averages. Slope ≈ 5.0; sub-5th order is recovered in the
  asymptotic regime.</figcaption>
</figure>

The numeric coverage lives in the Layer-B convergence fixture —
[<code>discretizations/finite_volume/weno5_advection_2d/fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/weno5_advection_2d/fixtures/convergence)
— which declares <code>expected_min_order = 4.5</code> on
<code>n ∈ {32, 64, 128, 256}</code> with the separable phase-shifted sine
product MMS <code>u(x,y) = sin(2π x + 1)·sin(2π y + 1)</code> on
<code>[0,1]²</code> periodic, matching the floor cited by Shu (1998) §2.2
with a small allowance for the ε-regularised nonlinear-weight transition.

Theoretical asymptotic order: **5.0** per axis (Jiang & Shu 1996, smooth
regions; Shu 1998 §2.2 dimension-by-dimension splitting; LeVeque 2002 §20).
Acceptance threshold: **min(observed order) ≥ 4.5** across the
<code>n ∈ {32, 64, 128, 256}</code> sweep.

The sweep is dispatched through ESS's
<code>mms_weno5_convergence</code> 2D path (esm-hsa): the harness applies
the rule's <code>axes.x</code> 1D sub-stencil row-wise on the cell-averaged
manufactured field and compares the reconstructed face values against the
perpendicular-axis-averaged analytic face truth. WENO5's positive
homogeneity makes the result symmetric in <code>x ↔ y</code>, so verifying
one axis is sufficient for the asymptotic-order claim.
