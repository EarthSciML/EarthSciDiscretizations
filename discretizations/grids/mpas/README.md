# MPAS grid fixtures

Canonical declarative `.esm` grid instances for the MPAS (unstructured Voronoi)
family. Fixtures live beside the family schema
(`discretizations/grids/mpas.schema.json`) and name standard resolutions from
the MPAS atmosphere mesh library.

Per the 2026-04-20 mayor scope correction on bead `dsc-3nw`, a `.esm` grid
entry is a **small declarative config** (family + options + dimensions +
loader ref), not a serialized geometry blob. The per-binding runtime provides
accessors (`cell_centers(c)`, `neighbors(c)`, `metric_eval(name, i)`,
`cell_area(c)`, `edge_length(e)`) that derive geometry on demand from the
referenced mesh file. Cross-binding conformance is verified by comparing
accessor outputs at pinned query points per `docs/GRIDS_API.md` §4.

## Fixtures

| File              | Mesh ID    | `n_cells` | `n_edges` | `n_vertices` | Nominal spacing |
|-------------------|------------|-----------|-----------|--------------|-----------------|
| `x1.642.esm`      | `x1.642`   | 642       | 1920      | 1280         | ~480 km         |
| `x1.2562.esm`     | `x1.2562`  | 2562      | 7680      | 5120         | ~240 km         |
| `x1.10242.esm`    | `x1.10242` | 10242     | 30720     | 20480        | ~120 km         |

All fixtures are quasi-uniform icosahedral-dual global Voronoi meshes (12
pentagons + `n_cells - 12` hexagons, `max_edges = 6`). Dimensions obey the
closed-sphere Euler identity for an icosahedral triangulation:

```
n_edges    = 3 * (n_cells - 2)
n_vertices = 2 * (n_cells - 2)
```

## Shape

Each `.esm` file has this structure (matches the binding `to_esm()` /
`toESM()` lowering minus the binding-specific `provenance.binding` tag —
fixtures are binding-neutral):

```jsonc
{
  "family": "mpas",
  "version": "1.0.0",
  "schema_version": "1.0.0",
  "topology": "unstructured",
  "options": {
    "R": 6371000.0,
    "dtype": "float64",
    "ghosts": 0,
    "loader": {
      "path": "meshes/mpas/<mesh_id>.grid.nc",
      "reader": "mpas_mesh",
      "check": "strict"
    }
  },
  "dimensions": {
    "n_cells": <int>,
    "n_edges": <int>,
    "n_vertices": <int>,
    "max_edges": 6
  },
  "provenance": {
    "source": "<human-readable mesh origin>",
    "mesh_id": "<x1.N>",
    "family": "mpas",
    "version": "1.0.0",
    "notes": "<string>"
  }
}
```

## Loader paths

`options.loader.path` is a repo-relative convention (`meshes/mpas/<mesh>.grid.nc`);
NetCDF mesh files are not bundled with this repo. Consumers supply a
`reader_fn` at bind time that translates the declared path into an in-memory
mesh per `docs/GRIDS_API.md` §10 — see, e.g., the Python
`earthsci_toolkit.grids.mpas(..., reader_fn=...)` signature and the TypeScript
`readerFn` contract in `typescript/src/grids/mpas.ts`.

## Regeneration

The committed dimensions are derived from the MPAS atmosphere mesh library
(`https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html`). If the mesh
library revises a standard resolution's cell count, regenerate by reading the
source mesh and updating the corresponding `.esm` dimensions block; then
re-run the Julia walker (`test/test_mpas_fixtures.jl`) to confirm the closed-
sphere Euler relations still hold.
