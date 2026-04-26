---
title: "Vertical"
slug: "vertical"
grid_families: "vertical"
rule_kinds: "grid"
description: "1D vertical column with arbitrary level spacing — pressure / sigma / hybrid coordinates."
source: "src/grids/vertical.jl"
tags: ["grid", "structured", "vertical", "column", "sigma", "hybrid"]
---

## Description

The `vertical` family is a stand-alone 1D column of stacked levels. It is the
column factor in any 2D-horizontal × 1D-vertical product grid, and it is also
used directly by box and column models.

Level edges may be supplied in any of the standard atmospheric coordinates:

- **z** (altitude in metres)
- **p** (pressure in Pa)
- **σ = (p − pₜ) / (pₛ − pₜ)** (terrain-following)
- **hybrid η = A(η) + B(η) (pₛ − pₜ)** (CAM/IFS-style)

The accessor stores the per-level edge array directly; cell widths
(`Δz`, `Δp`, `Δσ`, …) are derived. The trait registration is against
`AbstractVerticalGrid`.

## Visualization

<figure class="figure">
  <img src="/plots/grids/vertical.png" alt="Stretched vertical levels">
  <figcaption>12 stretched vertical levels — sigma-like spacing, dense near
  the surface and sparse aloft.</figcaption>
</figure>

## Trait coverage

`AbstractVerticalGrid`. Implements `cell_centers`, `cell_widths`,
`neighbor_indices` over the single vertical axis.

## See also

- [`centered_2nd_uniform_vertical`]({{< ref "/rules/centered_2nd_uniform_vertical" >}})
- Lin (2004), MWR — vertical Lagrangian remap
- [GRIDS_API.md]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md)
