---
title: "divergence_arakawa_c"
slug: "divergence_arakawa_c"
family: "finite_volume"
grid_family: "arakawa"
rule_kind: "scheme"
accuracy: "O(h²)"
applies_to: "div(F)"
rule_path: "discretizations/finite_volume/divergence_arakawa_c.json"
description: "Two-point centered finite-volume divergence on a C-grid — F_x at face_x, F_y at face_y, divergence emitted at cell_center."
tags: ["finite-volume", "arakawa", "c-grid", "divergence", "staggered"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/divergence_arakawa_c-stencil.png"
       alt="C-grid divergence — 4-point stencil over the cell faces">
  <figcaption>The C-grid divergence reads <code>F<sub>x</sub></code> at the two
  face_x edges and <code>F<sub>y</sub></code> at the two face_y edges, and emits
  <code>div F</code> at the cell center. Stagger metadata is carried in the rule
  via <code>requires_locations</code> and <code>emits_location</code>.</figcaption>
</figure>

## Coefficients

| selector | stagger | axis | offset | coeff |
|---|---|---|---:|---|
| `arakawa` | `face_x` | `$x` |  0 | `−1 / dx` |
| `arakawa` | `face_x` | `$x` | +1 | `+1 / dx` |
| `arakawa` | `face_y` | `$y` |  0 | `−1 / dy` |
| `arakawa` | `face_y` | `$y` | +1 | `+1 / dy` |

The rule declares an `arakawa_stagger` enum mapping stagger names
(`cell_center`, `face_x`, `face_y`, `vertex`) to integer codes — the
runtime uses these to dispatch `location_centers(grid, ...)` per stagger.

`requires_locations: ["face_x", "face_y"]`,
`emits_location: "cell_center"`.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B harness needs stagger-aware sampling — manufactured solutions
must be evaluated at the appropriate <code>face_x</code> / <code>face_y</code> location for each
input field, not just at cell centers. The fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/divergence_arakawa_c/fixtures/convergence)
declares <code>applicable: false</code> until that lands.
</div>
