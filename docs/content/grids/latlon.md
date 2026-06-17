---
title: "Lat-Lon"
slug: "latlon"
grid_families: "latlon"
rule_kinds: "grid"
description: "Regular spherical lat-lon mesh with metric corrections for the spherical surface."
source: "src/grids/latlon.jl"
tags: ["grid", "structured", "spherical", "lat-lon"]
---

## Description

The lat-lon family is a regular mesh on the sphere of radius `R` indexed by
longitude λ and latitude φ (matching the code's metric convention in
`src/grids/latlon.jl`). Cells are bounded by meridians (constant λ) and
parallels (constant φ); cell area shrinks toward the poles as `R² cos φ dφ dλ`.
The metric tensor is diagonal with `g_{λλ} = R² cos² φ` and `g_{φφ} = R²`,
so finite-difference rules carry a `cos φ` factor on the longitudinal stencil
(see [`centered_2nd_uniform_latlon`]({{< ref "/rules/centered_2nd_uniform_latlon" >}})).

Pole handling is configurable via the `pole_policy` knob — `none` (the default)
leaves polar cells degenerate, while wrap-and-stitch policies live in the
accessor.

## Visualization

<figure class="figure">
  <img src="/plots/grids/latlon.png" alt="Lat-lon mesh in equirectangular projection">
  <figcaption>Lat-lon 24×12 mesh shown in equirectangular projection.
  Cell area is proportional to cos φ (latitude) — visible in the spherical
  metric, not in this flat projection.</figcaption>
</figure>

## Trait coverage

Registered against `AbstractCurvilinearGrid`. Supplies the full bulk-array
trait contract; `metric_g` / `metric_ginv` / `metric_jacobian` carry the
spherical metric, which downstream rules consume directly.

## Canonical fixtures

- `discretizations/grids/latlon/lat_lon_1deg.esm` — 360×180 (1°)
- `discretizations/grids/latlon/lat_lon_0p25deg.esm` — 1440×720 (0.25°)
- `discretizations/grids/latlon/lat_lon_0p1deg.esm` — 3600×1800 (0.1°)

## See also

- [`centered_2nd_uniform_latlon`]({{< ref "/rules/centered_2nd_uniform_latlon" >}})
- [GRIDS_API.md]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md)
