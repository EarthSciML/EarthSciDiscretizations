---
title: "flux_limiter_minmod"
slug: "flux_limiter_minmod"
family: "limiter"
grid_family: "cartesian"
rule_kind: "limiter"
accuracy: "O(dx²) in smooth monotone regions; O(dx) at extrema"
applies_to: "limit(r)"
rule_path: "discretizations/finite_volume/flux_limiter_minmod.json"
description: "Minmod TVD flux limiter (Roe 1986) — the most diffusive of the symmetric second-order TVD limiters."
tags: ["limiter", "tvd", "minmod", "monotonicity"]
---

## Limiter curve

<figure class="figure">
  <img src="/plots/rules/flux_limiter_minmod-stencil.png"
       alt="Minmod limiter φ(r) and Sweby second-order TVD region">
  <figcaption>Minmod limiter <code>φ(r) = max(0, min(r, 1))</code> overlaid on the
  Sweby second-order TVD region. The curve sits on the lower edge of the
  region and passes through <code>(1, 1)</code>.</figcaption>
</figure>

## Formula

$$\varphi(r) = \max(0, \min(r, 1))$$

with `r = (q_i − q_{i−1}) / (q_{i+1} − q_i)` evaluated at each interface.
Properties:

| Property | Value |
|---|---|
| TVD | yes |
| Monotonicity-preserving | yes |
| Symmetric | yes |
| Smooth-extremum behavior | drops to 1st-order by design |

Caller multiplies the high-order slope correction by `φ(r)` at each edge —
see the worked example in
[`discretizations/finite_volume/README.md`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/README.md).

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
Layer-B' (monotonicity / TVD) is the right harness kind for limiters — Layer-B
convergence does not apply directly because limiters intentionally degrade to
first order at smooth extrema. The
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_minmod/fixtures/convergence)
fixture is provisional; once the limiter MMS harness is wired up (see
walker rewrite/ and limiter/ fixture kinds), this section will show
TVD-norm tests instead of order-of-convergence plots.
</div>
