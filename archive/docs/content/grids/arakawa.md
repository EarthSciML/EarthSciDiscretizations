---
title: "Arakawa staggered grid"
slug: "arakawa"
grid_families: "arakawa"
rule_kinds: "grid"
description: "Staggered grid wrapping a base curvilinear grid — A / B / C / D / E variants — with per-stagger variable locations."
source: "src/grids/arakawa.jl"
tags: ["grid", "staggered", "arakawa", "c-grid", "b-grid", "shallow-water"]
---

## Description

The Arakawa family is *not* a standalone topology — it is a stagger wrapper
around a base curvilinear grid (Cartesian or lat-lon). It assigns prognostic
variables to one of four canonical
locations per cell. The location enum (`VarLocation`, exported from the
package) is `CellCenter`, `UEdge`, `VEdge`, `Corner`:

| `VarLocation` | Accessor | Variable held there | Used by |
|---|---|---|---|
| `CellCenter` | `cell_centers(grid, i, j)` | scalars (mass `h`, pressure, tracers) | A, B, C, D, E |
| `UEdge` | `u_face(grid, i, j)` | u-component of velocity | C, D, E |
| `VEdge` | `v_face(grid, i, j)` | v-component of velocity | C, D, E |
| `Corner` | `corners(grid, i, j)` | u, v together (B) | B |

The Arakawa A/B/C/D/E distinction is just the choice of which subset of
locations carries `u`, `v`, and `h` (the per-stagger table lives in
`arakawa_variable_locations` in `src/grids/arakawa.jl`):

- **A-grid** (`ArakawaA`): all variables co-located at cell centers (simple
  but hosts checkerboard modes for the divergence operator)
- **B-grid** (`ArakawaB`): `h` at center, both `u` and `v` at corners
- **C-grid** (`ArakawaC`): `h` at center, `u` on the U-edge, `v` on the
  V-edge (good for wave equations and divergence-conservative
  discretizations)
- **D-grid** (`ArakawaD`): dual of C — `u` and `v` swapped between the
  U-edge and V-edge
- **E-grid** (`ArakawaE`): rotated arrangement with `u` and `v` on edges
  (the `stagger` parser accepts `:E`)

A grid is constructed via
`EarthSciDiscretizations.grids.arakawa(; base, stagger, ghosts, dtype)`,
where `stagger` is one of `:A`, `:B`, `:C`, `:D`, `:E` and `base` is the
underlying Cartesian or lat-lon grid.

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
underlying Cartesian / lat-lon implementation; stagger-aware accessors —
`cell_centers`, `u_face`, `v_face`, `corners` — return the coordinate at
the location where each prognostic variable lives for the chosen stagger.

## See also

- [`divergence_arakawa_c`]({{< ref "/rules/divergence_arakawa_c" >}})
- Arakawa & Lamb (1977), Methods in Computational Physics 17:173–265
- `discretizations/SELECTOR_KINDS.md`
