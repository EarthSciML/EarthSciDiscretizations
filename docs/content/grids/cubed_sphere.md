---
title: "Cubed-sphere"
slug: "cubed_sphere"
grid_families: "cubed_sphere"
rule_kinds: "grid"
description: "Six gnomonic panels covering the sphere; in-panel curvilinear (ξ, η) coordinates with cross-panel ghost stitching."
source: "src/grids/cubed_sphere.jl"
tags: ["grid", "block-structured", "spherical", "cubed-sphere", "gnomonic"]
---

## Description

A cubed-sphere mesh tiles the sphere with six logically rectangular panels —
one per face of an inscribed cube — projected outward by a gnomonic
projection. Each panel carries an in-panel `(ξ, η)` curvilinear coordinate
system on `[-1, 1]²`, discretized into `Nc × Nc` cells with uniform isotropic
spacing `h = π/(2 Nc)`.

The metric tensor is non-diagonal: `g^{ξξ}`, `g^{ηη}`, `g^{ξη}` are all non-zero
and depend on position within the panel. Cross-panel adjacency is resolved by
the panel-connectivity table (`src/grids/panel_connectivity.jl`); discretization
rules reference only in-panel offsets and rely on the accessor to fill cross-panel
ghosts and rotate vector components on panel boundaries.

## Visualization

<figure class="figure">
  <img src="/plots/grids/cubed_sphere.png" alt="Cubed-sphere C6 mesh">
  <figcaption>Gnomonic cubed-sphere with `Nc = 6` (216 cells across 6 panels).
  Each panel is colored differently. The non-uniform pixel density visible at
  panel corners reflects the well-known gnomonic distortion.</figcaption>
</figure>

## Trait coverage

Subtype of `AbstractCubedSphereGrid` (which is itself a curvilinear grid).
Implements the bulk-array contract panel-by-panel; `metric_g`, `metric_ginv`,
`metric_jacobian`, and `metric_dgij_dxk` are required by the
[`covariant_laplacian_cubed_sphere`]({{< ref "/rules/covariant_laplacian_cubed_sphere" >}})
rule.

## Canonical fixtures

- `discretizations/grids/cubed_sphere/c24.esm` — Nc=24
- `discretizations/grids/cubed_sphere/c48.esm` — Nc=48
- `discretizations/grids/cubed_sphere/c96.esm` — Nc=96

## See also

- [`covariant_laplacian_cubed_sphere`]({{< ref "/rules/covariant_laplacian_cubed_sphere" >}})
- Putman & Lin (2007), JCP 227(1):55–78 — gnomonic cubed-sphere covariant FV operators
