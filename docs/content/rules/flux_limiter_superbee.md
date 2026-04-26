---
title: "flux_limiter_superbee"
slug: "flux_limiter_superbee"
family: "limiter"
grid_family: "cartesian"
rule_kind: "limiter"
accuracy: "O(dx²) in smooth monotone regions; compressive near discontinuities"
applies_to: "limit(r)"
rule_path: "discretizations/finite_volume/flux_limiter_superbee.json"
description: "Superbee TVD flux limiter (Roe 1986) — sits on the upper edge of the second-order TVD region."
tags: ["limiter", "tvd", "superbee", "compressive"]
---

## Limiter curve

<figure class="figure">
  <img src="/plots/rules/flux_limiter_superbee-stencil.png"
       alt="Superbee limiter φ(r) and Sweby second-order TVD region">
  <figcaption>Superbee limiter on the upper edge of the Sweby second-order
  TVD region — compressive: it amplifies sharp gradients but degrades
  smooth extrema more aggressively than minmod.</figcaption>
</figure>

## Formula

$$\varphi(r) = \max\bigl(0,\;\min(2r, 1),\;\min(r, 2)\bigr)$$

| Property | Value |
|---|---|
| TVD | yes |
| Monotonicity-preserving | yes |
| Symmetric | yes |
| Smooth-extremum behavior | degrades more than minmod by design |

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
As with [<code>flux_limiter_minmod</code>]({{< ref "/rules/flux_limiter_minmod" >}}),
the right Layer-B' (monotonicity / TVD) harness kind is in flight. The
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_superbee/fixtures/convergence)
fixture is provisional; final form will exercise TVD norms across initial
conditions with sharp gradients, not point-wise convergence.
</div>
