---
title: "MPAS"
slug: "mpas"
grid_families: "mpas"
rule_kinds: "grid"
description: "Quasi-uniform Voronoi mesh on the sphere with icosahedral dual (12 pentagons + hexagons)."
source: "src/grids/mpas.jl"
tags: ["grid", "unstructured", "voronoi", "spherical", "mpas", "soca"]
---

## Description

MPAS (Model for Prediction Across Scales) uses a spherical centroidal
Voronoi tessellation: a quasi-uniform mesh of polygonal cells (almost all
hexagons, with exactly 12 pentagons inherited from the icosahedral seed).
Cells are described not by `(i, j)` block-structured indices but by edge-
and vertex-of-cell adjacency tables loaded from a NetCDF mesh file:

| Table | Meaning |
|---|---|
| `cells_on_cell[c, k]` | k-th neighbor of cell `c` |
| `edges_on_cell[c, k]` | k-th edge of cell `c` |
| `n_edges_on_cell[c]`  | number of edges (5 or 6) |
| `dv_edge[e]`          | edge length (distance between Voronoi vertices) |
| `dc_edge[e]`          | dual-edge length (distance between cell centers) |
| `area_cell[c]`        | cell area |

Discretization rules on this family use indirect indexing (`reduction` and
`indirect` selectors — see `discretizations/SELECTOR_KINDS.md`) since
neighbor count varies per cell.

## Visualization

<figure class="figure">
  <img src="/plots/grids/mpas.png" alt="Quasi-uniform spherical Voronoi mesh">
  <figcaption>Illustrative spherical Voronoi mesh (162 cells). Real MPAS
  meshes (e.g. <code>x1.642</code>, <code>x1.10242</code>, <code>x1.40962</code>) are derived from
  successively-bisected icosahedral seeds and have the characteristic
  12 pentagons + (n − 12) hexagons.</figcaption>
</figure>

## Trait coverage

Registered against `AbstractUnstructuredGrid`. The `neighbor_indices` API
returns ragged adjacency arrays; bulk fields like `area_cell` and the per-edge
geometry are loaded once from the source NetCDF and exposed as flat arrays
keyed by cell / edge / vertex index.

## Canonical fixtures

- `discretizations/grids/mpas/x1.642.esm` — 642 cells (~480 km)
- `discretizations/grids/mpas/x1.2562.esm` — 2562 cells
- `discretizations/grids/mpas/x1.10242.esm` — 10242 cells (~120 km)

Mesh files themselves live outside the repo (the `.esm` is a manifest;
the loader resolves `meshes/mpas/*.grid.nc` at bind time per `GRIDS_API.md` §10).

## See also

- [`nn_diffusion_mpas`]({{< ref "/rules/nn_diffusion_mpas" >}})
- Skamarock et al. (2012), MWR 140
- Thuburn et al. (2009), JCP 228; Ringler et al. (2010), JCP 229
