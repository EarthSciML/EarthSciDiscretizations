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

The Layer-B walker drives this rule through ESS's
<code>verify_mms_convergence</code> with the registered
<code>Y<sub>2,0</sub></code> spherical-harmonic manufactured solution on the
unit sphere (lon-independent, so the lon stencil is exact and the test
isolates the lat-axis 2nd-order accuracy). The fixture sweeps
<code>n ∈ {16, 32, 64, 128}</code> with <code>nlon = 2 n</code>,
<code>nlat = n</code> and asserts an observed L<sub>∞</sub> minimum order
≥ 1.9 over interior cells (poles excluded by the lat stencil's reach).
Inputs/expected live under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_uniform_latlon/fixtures/convergence).
