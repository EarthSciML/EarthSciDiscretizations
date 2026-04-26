---
title: "Arakawa staggered grid"
slug: "arakawa"
grid_families: "arakawa"
rule_kinds: "grid"
description: "Staggered grid wrapping a base curvilinear grid — A / B / C / D variants — with per-stagger location keys."
source: "src/grids/arakawa.jl"
tags: ["grid", "staggered", "arakawa", "c-grid", "b-grid", "shallow-water"]
---

## Description

The Arakawa family is *not* a standalone topology — it is a stagger wrapper
around a base curvilinear grid (Cartesian or lat-lon today; cubed-sphere is
on the roadmap). It assigns prognostic variables to one of four canonical
locations per cell:

| Location | Variable held there | Used by |
|---|---|---|
| `cell_center` | scalars (mass `h`, pressure, tracers) | A, B, C, D |
| `face_x` | u-component of velocity (C/D) | C |
| `face_y` | v-component of velocity (C/D) | C |
| `vertex` | u, v together (B) | B |

The Arakawa A/B/C/D distinction is just the choice of which subset of
locations carries `u`, `v`, and `h`:

- **A-grid**: all variables co-located at cell centers (simple but
  hosts checkerboard modes for the divergence operator)
- **B-grid**: `h` at center, both `u` and `v` at vertices
- **C-grid**: `h` at center, `u` at face_x, `v` at face_y (good for
  wave equations and divergence-conservative discretizations)
- **D-grid**: dual of C — `u` and `v` swapped between face_x and face_y

ESD currently seeds the C-grid via the
[`divergence_arakawa_c`]({{< ref "/rules/divergence_arakawa_c" >}}) rule, which
demonstrates the `stagger_enum` / `requires_locations` / `emits_location`
fields used to make stagger-position selectors authoritative in the rule
schema.

## Visualization

<figure class="figure">
  <img src="/plots/grids/arakawa.png" alt="Arakawa A/B/C/D stagger comparison">
  <figcaption>Variable placement on a single cell for each Arakawa stagger:
  A (all co-located), B (uv on vertices), C (u on face_x, v on face_y),
  D (the dual of C).</figcaption>
</figure>

## Trait coverage

Registered against `AbstractStaggeredGrid` (which sits beneath
`AbstractCurvilinearGrid`). The base-grid trait calls all delegate to the
underlying Cartesian / lat-lon implementation; stagger-aware methods like
`location_centers(grid, :face_x)` add the per-location offset.

## See also

- [`divergence_arakawa_c`]({{< ref "/rules/divergence_arakawa_c" >}})
- Arakawa & Lamb (1977), Methods in Computational Physics 17:173–265
- `discretizations/SELECTOR_KINDS.md`
