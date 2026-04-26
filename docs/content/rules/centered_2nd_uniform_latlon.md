---
title: "centered_2nd_uniform_latlon"
slug: "centered_2nd_uniform_latlon"
family: "finite_difference"
grid_family: "latlon"
rule_kind: "scheme"
accuracy: "O(h²)"
applies_to: "grad(u), dim=k (k ∈ {lon, lat})"
rule_path: "discretizations/finite_difference/centered_2nd_uniform_latlon.json"
description: "Centered second-order finite difference on a regular lat-lon grid, with the standard cos λ metric on the longitudinal stencil."
tags: ["finite-difference", "lat-lon", "spherical", "metric"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/centered_2nd_uniform_latlon-stencil.png"
       alt="centered_2nd_uniform_latlon stencils for the lon and lat axes">
  <figcaption>Two-point centered stencil applied independently on the lon and
  lat axes; the lon coefficients carry the <code>cos λ</code> spherical metric
  factor.</figcaption>
</figure>

## Coefficients

| selector kind | axis  | offset | coeff |
|---|---|---:|---|
| `latlon` | `lon` | −1 | `−1 / (2 R cos λ dλ)` |
| `latlon` | `lon` | +1 | `+1 / (2 R cos λ dλ)` |
| `latlon` | `lat` | −1 | `−1 / (2 R dφ)` |
| `latlon` | `lat` | +1 | `+1 / (2 R dφ)` |

`R` is sphere radius, `λ` is latitude, `dλ` and `dφ` are constant grid
spacings. The lat-lon accessor surfaces `cos_lat` as a per-row metric
quantity bound at evaluation time.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B convergence harness needs the in-flight 2D dispatch and per-cell
metric callables (see
<code>EarthSciSerialization/packages/EarthSciSerialization.jl/src/mms_evaluator.jl</code>)
to evaluate this rule on a spherical manufactured solution. The fixture
under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_latlon/fixtures/convergence)
declares <code>applicable: false</code> until those extensions land.
</div>
