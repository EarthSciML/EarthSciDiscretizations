---
title: "covariant_laplacian_cubed_sphere"
slug: "covariant_laplacian_cubed_sphere"
families: "finite_difference"
grid_families: "cubed_sphere"
rule_kinds: "scheme"
accuracy: "O(h²)"
applies_to: "laplacian(φ)"
rule_path: "discretizations/finite_difference/covariant_laplacian_cubed_sphere.json"
description: "9-point covariant Laplacian on the gnomonic cubed sphere — expressed as a cross_metric composite (§7.6) with 8 terms covering orthogonal, cross-derivative, and metric-gradient correction contributions."
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

## Continuous form

Starting from the divergence-form Laplace–Beltrami operator on the gnomonic
cubed-sphere panel,

$$\nabla^2 \varphi = \tfrac{1}{J}\bigl[\partial_\xi(J\,g^{\xi\xi}\,\partial_\xi\varphi + J\,g^{\xi\eta}\,\partial_\eta\varphi) + \partial_\eta(J\,g^{\xi\eta}\,\partial_\xi\varphi + J\,g^{\eta\eta}\,\partial_\eta\varphi)\bigr],$$

we expand into an **orthogonal part** plus a **cross-metric part**:

$$\nabla^2 \varphi = \underbrace{g^{\xi\xi}\,\partial^2_\xi\varphi + g^{\eta\eta}\,\partial^2_\eta\varphi + \tfrac{1}{J}\bigl(\partial_\xi(J g^{\xi\xi})\,\partial_\xi\varphi + \partial_\eta(J g^{\eta\eta})\,\partial_\eta\varphi\bigr)}_{\text{orthogonal}} + \underbrace{2\,g^{\xi\eta}\,\partial^2_{\xi\eta}\varphi + \tfrac{1}{J}\bigl(\partial_\xi(J g^{\xi\eta})\,\partial_\eta\varphi + \partial_\eta(J g^{\xi\eta})\,\partial_\xi\varphi\bigr)}_{\text{cross-metric}}.$$

## Discrete operator

Substituting centered second-order finite differences for the first and
second derivatives over in-panel cell offsets `(Δξ, Δη) ∈ {−1, 0, +1}²`
yields a 9-point linear combination of cell values:

$$(\nabla^2\varphi)_{i,j} \;=\; \sum_{a,b\in\{-1,0,+1\}} c_{a,b}\;\varphi_{i+a,\,j+b},$$

with coefficients (uniform isotropic spacing `h = dξ = dη = π/(2 N_c)`):

| offset `(Δξ, Δη)` | coefficient `c_{Δξ,Δη}` | contribution |
|---|---|---|
| `( 0,  0)` | `−2 (g^{ξξ} + g^{ηη}) / h²` | central — orthogonal `∂²_ξ + ∂²_η` diagonal |
| `(+1,  0)` | `g^{ξξ}/h² + (∂_ξ(J g^{ξξ}) + ∂_η(J g^{ξη})) / (2 h J)` | east — orthogonal `∂²_ξ` + first-derivative correction |
| `(−1,  0)` | `g^{ξξ}/h² − (∂_ξ(J g^{ξξ}) + ∂_η(J g^{ξη})) / (2 h J)` | west |
| `( 0, +1)` | `g^{ηη}/h² + (∂_η(J g^{ηη}) + ∂_ξ(J g^{ξη})) / (2 h J)` | north — orthogonal `∂²_η` + first-derivative correction |
| `( 0, −1)` | `g^{ηη}/h² − (∂_η(J g^{ηη}) + ∂_ξ(J g^{ξη})) / (2 h J)` | south |
| `(+1, +1)` | `+g^{ξη} / (2 h²)` | NE corner — cross-metric `∂²_{ξη}` |
| `(−1, −1)` | `+g^{ξη} / (2 h²)` | SW corner |
| `(+1, −1)` | `−g^{ξη} / (2 h²)` | SE corner |
| `(−1, +1)` | `−g^{ξη} / (2 h²)` | NW corner |

Each metric quantity is resolved at the **central cell** `(i, j)` by the
`cubed_sphere` grid accessor; the rule does not interpolate metric to face
midpoints. The full machine-readable coefficient expressions live in
[the JSON rule file]({{< param repoURL >}}/blob/main/discretizations/finite_difference/covariant_laplacian_cubed_sphere.json).

### Metric bindings

| Binding | Resolved by |
|---|---|
| `J`, `ginv_xi_xi`, `ginv_eta_eta`, `ginv_xi_eta` | cubed_sphere accessor at the central cell |
| `dJgxx_dxi`, `dJgyy_deta`, `dJgxe_dxi`, `dJgxe_deta` | metric-derivative tables on the panel |
| `h` | uniform isotropic spacing `dξ = dη = π / (2 Nc)` |

Cross-panel ghost cells and basis rotation at panel boundaries live in the
cubed_sphere accessor (`src/grids/panel_connectivity.jl`); the per-axis
stencils within the cross_metric composite do **not** carry a `panel` field —
`(Δξ, Δη)` offsets resolve across panel seams via the accessor's
connectivity table.

## Convergence

Numeric coverage today lives in the canonical Julia test
the walker's `cubed_sphere_cross_metric` Layer-B runner (esd-ecq):
the manufactured solution `φ(ξ, η) = cos(2ξ)·cos(2η)` is evaluated against the
analytic covariant Laplacian (with metric quantities sampled from the gnomonic
metric) at `Nc ∈ {64, 128, 256}`. Coarser grids (Nc < 32) show pre-asymptotic
convergence (~O(h^1.2)) due to the gnomonic metric's large variation near panel
corners. The asymptotic O(h²) regime is reached at Nc ≥ 64 (order ≥ 1.9).
The convergence fixture at
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/covariant_laplacian_cubed_sphere/fixtures/convergence)
declares `mms_kind="cos2xi_cos2eta_cubed_sphere"` and asserts `expected_min_order: 1.8`.
The imperative oracle `fv_laplacian` (`src/operators/laplacian.jl`) was retired
in esd-ecq; the declarative rule is the sole implementation.

The bit-equivalence canonical fixture at
[<code>fixtures/canonical/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/covariant_laplacian_cubed_sphere/fixtures/canonical)
remains `applicable: false` pending ESS cubed_sphere selector dispatch (tracked
at ESS/ess-cubed-sphere-pipeline).

## Reference

- Putman & Lin (2007), *JCP* 227(1):55–78 — gnomonic cubed-sphere covariant FV operators.
