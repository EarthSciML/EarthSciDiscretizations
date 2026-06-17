---
title: "Vertical"
slug: "vertical"
grid_families: "vertical"
rule_kinds: "grid"
description: "1D vertical column with arbitrary level spacing — sigma / eta / z / theta / hybrid / z_star coordinates."
source: "src/grids/vertical.jl"
tags: ["grid", "structured", "vertical", "column", "sigma", "hybrid"]
---

## Description

The `vertical` family is a stand-alone 1D column of stacked levels. It is the
column factor in any 2D-horizontal × 1D-vertical product grid, and it is also
used directly by box and column models.

The `coordinate` option selects the vertical coordinate. The accepted kinds
(see `_VERTICAL_COORDINATES` in `src/grids/vertical.jl`) are:

- **`:z`** — altitude in metres
- **`:sigma`** — σ = (p − pₜ) / (pₛ − pₜ), terrain-following
- **`:eta`** — η hybrid-coordinate parameter
- **`:theta`** — potential-temperature (isentropic) coordinate
- **`:hybrid_sigma_theta`** — hybrid σ–θ coordinate (`ak`/`bk` coefficients)
- **`:z_star`** — height-based terrain-following (z\*) coordinate

(There is no `:pressure` coordinate — `pressure` is a *metric name*, derived
from the levels, not a `coordinate` value.) The accessor stores the per-level
edge array directly; cell widths are derived. The trait registration is
against `AbstractVerticalGrid`.

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
