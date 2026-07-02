---
title: "Duo (geodesic triangular)"
slug: "duo"
grid_families: "duo"
rule_kinds: "grid"
description: "Geodesic triangular mesh on the sphere; recursively subdivided icosahedron."
source: "src/grids/duo.jl"
tags: ["grid", "unstructured", "triangular", "icosahedral", "geodesic"]
---

## Description

The `duo` family is a geodesic triangular mesh — a recursively-bisected
icosahedron projected onto the sphere. At level `L`, each of the 20
icosahedral seed faces is subdivided into 4ᴸ triangles, giving `20·4ᴸ`
triangles, `30·4ᴸ` edges, and `2 + 10·4ᴸ` vertices.

The mesh is unstructured (no `(i, j)` indexing) but enjoys topological
regularity: every interior vertex has 6 incident edges and the only
"defects" are the 12 vertices inherited from the icosahedron, which have
5 edges. This makes it a popular target for shallow-water and barotropic
codes that want a uniform sphere coverage without the cubed-sphere panel
seams.

## Visualization

<figure class="figure">
  <img src="/plots/grids/duo.png" alt="Geodesic triangular mesh, level 1">
  <figcaption>Level-1 icosahedral subdivision (80 triangles, 42 vertices).
  The 12 original icosahedron vertices remain visible as the only
  five-coordinated nodes.</figcaption>
</figure>

## Trait coverage

Registered against `AbstractUnstructuredGrid`. Adjacency arrays are
`vertices_on_face`, `edges_on_face`, and `faces_on_vertex`; bulk geometry
is the per-vertex Cartesian-on-sphere position and per-edge arc length.

## Canonical fixtures

- `discretizations/grids/duo/icos_level0.esm` — 20 faces, 12 vertices (the seed)
- `discretizations/grids/duo/icos_level1.esm` — 80 faces, 42 vertices
- `discretizations/grids/duo/icos_level2.esm` — 320 faces, 162 vertices

## See also

- [GRIDS_API.md]({{< param repoURL >}}/blob/main/docs/GRIDS_API.md)
