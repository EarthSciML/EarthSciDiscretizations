---
title: "covariant_laplacian_cubed_sphere"
slug: "covariant_laplacian_cubed_sphere"
family: "finite_difference"
grid_family: "cubed_sphere"
rule_kind: "scheme"
accuracy: "O(h²)"
applies_to: "laplacian(φ)"
rule_path: "discretizations/finite_difference/covariant_laplacian_cubed_sphere.json"
description: "9-point covariant Laplacian on the gnomonic cubed sphere — orthogonal + cross-metric corrections from the inverse metric tensor."
tags: ["finite-difference", "cubed-sphere", "covariant", "laplacian"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/covariant_laplacian_cubed_sphere-stencil.png"
       alt="9-point in-panel stencil with metric-weighted coefficients">
  <figcaption>9-point stencil over the in-panel <code>(ξ, η)</code> coordinates. The
  4-point cross stencil carries the orthogonal metric (<code>g^{ξξ}</code>, <code>g^{ηη}</code>)
  and metric-derivative corrections; the 4 corner cells carry the
  cross-metric contribution <code>±g^{ξη} / (2 h²)</code>.</figcaption>
</figure>

## Coefficients

The full Laplacian is

$$\nabla^2 \varphi = \tfrac{1}{J}\bigl[\partial_\xi(J\,g^{\xi\xi}\,\partial_\xi\varphi + J\,g^{\xi\eta}\,\partial_\eta\varphi) + \partial_\eta(J\,g^{\xi\eta}\,\partial_\xi\varphi + J\,g^{\eta\eta}\,\partial_\eta\varphi)\bigr].$$

Expanded into orthogonal and cross-metric parts and discretized with a
9-point centered second-order stencil. The full coefficient table lives in
the rule file (see [the JSON]({{< param repoURL >}}/blob/main/discretizations/finite_difference/covariant_laplacian_cubed_sphere.json));
key bindings:

| Binding | Resolved by |
|---|---|
| `J`, `ginv_xi_xi`, `ginv_eta_eta`, `ginv_xi_eta` | cubed_sphere accessor at the central cell |
| `dJgxx_dxi`, `dJgyy_deta`, `dJgxe_dxi`, `dJgxe_deta` | metric-derivative tables on the panel |
| `h` | uniform isotropic spacing `dξ = dη = π / (2 Nc)` |

Cross-panel ghost cells and basis rotation at panel boundaries live in the
cubed_sphere accessor (`src/grids/panel_connectivity.jl`); the rule
selectors do **not** carry a `panel` field.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B harness needs the in-flight 2D dispatch + per-cell metric
callables to evaluate this rule on a manufactured solution defined in
<code>(ξ, η)</code> with the gnomonic Jacobian threaded through. The fixture
under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/covariant_laplacian_cubed_sphere/fixtures/convergence)
will populate once those extensions land.
</div>
