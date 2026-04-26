---
title: "Cartesian"
slug: "cartesian"
grid_family: "cartesian"
rule_kind: "grid"
description: "Logically rectangular Cartesian mesh with optional non-uniform spacing per axis (1D / 2D / 3D)."
source: "src/grids/cartesian.jl"
tags: ["grid", "structured", "rectilinear", "uniform", "non-uniform"]
---

## Description

The Cartesian grid is the simplest curvilinear topology ESD supports: a logically
rectangular mesh of cells whose coordinates are an outer product of per-axis edge
arrays. Each axis can be uniform (single-step `dx`) or non-uniform (arbitrary
1D edge array). Dimensionality is 1 / 2 / 3 (axes named `x`, `y`, `z`).

The on-disk `.esm` form is a small declarative config — `family`, `n[]`,
`extent[]`, `uniform[]` — not a serialized geometry blob. The runtime accessor
in `src/grids/cartesian.jl` derives `cell_centers`, `cell_widths`, `cell_volume`,
and `neighbor_indices` from that declaration via pure math.

## Visualization

<figure class="figure">
  <img src="/plots/grids/cartesian.png" alt="Cartesian 2D mesh">
  <figcaption>2D uniform Cartesian mesh, 8×6 cells. Yellow markers are
  cell centers; lines are cell edges.</figcaption>
</figure>

## Trait coverage

The Cartesian family registers against `AbstractCurvilinearGrid` and supplies
the full Tier-C bulk-array contract: `cell_centers`, `cell_widths`,
`cell_volume`, `neighbor_indices`, `boundary_mask`, `n_cells`, `n_dims`,
`axis_names`. Because the metric tensor is identity, `metric_g`, `metric_ginv`,
`metric_jacobian`, and `metric_dgij_dxk` short-circuit to constant arrays.

## Canonical fixtures

- `discretizations/grids/cartesian/uniform_1d_n64.esm` — N=64 uniform 1D mesh

## See also

- [`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}})
- [`upwind_1st`]({{< ref "/rules/upwind_1st" >}})
- [GRIDS_API.md]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md)
