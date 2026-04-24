# DUO grid fixtures

Canonical `.esm` declarative configurations for the `duo` icosahedral
triangular grid family (Heikes et al. 2023), per the cross-binding API
contract in `docs/GRIDS_API.md` §10.

Per the 2026-04-20 mayor correction on bead `dsc-8qn`, a `.esm` grid
entry is a **small declarative config** (family + options + loader ref +
dimensions), NOT a serialized geometry blob. Each binding's runtime
provides accessors (`cell_centers`, `neighbors`, `metric_eval`, `total_area`)
that DERIVE geometry from the declaration — for DUO, via the
`builtin://icosahedral/<level>` loader (pure math, no file I/O).

Cross-binding conformance is gated on accessor outputs at pinned cell
indices rather than serialized geometry bytes; see
`tests/conformance/grids/duo/` for that harness.

## Fixtures

Three canonical subdivision levels at Earth radius
(`R = 6.371e6 m`, `dtype = float64`, `ghosts = 0`). Subdivision level `r`
yields `20·4^r` triangular cells, `10·4^r + 2` vertices, `30·4^r` edges;
a closed icosahedral mesh satisfies the sphere Euler identity
`V - E + F = 2`.

| File                | level | `n_cells` | `n_vertices` | `n_edges` | Nominal spacing |
|---------------------|-------|-----------|--------------|-----------|-----------------|
| `icos_level0.esm`   |  0    |    20     |     12       |    30     | ~7000 km (bare icosahedron) |
| `icos_level1.esm`   |  1    |    80     |     42       |   120     | ~3500 km        |
| `icos_level2.esm`   |  2    |   320     |    162       |   480     | ~1800 km        |

`icos_level0` covers the bare icosahedron (exercises the generator base
case); `icos_level1` and `icos_level2` exercise recursive subdivision and
the midpoint cache.

## Shape

Each `.esm` file is binding-neutral and matches the Julia binding's
`to_esm(DuoGrid)` lowering minus the binding-specific `provenance` block
(per `docs/GRIDS_API.md` §6.4):

```jsonc
{
  "family": "duo",
  "topology": "unstructured",
  "dtype": "float64",
  "ghosts": 0,
  "n_cells": <int>,
  "n_vertices": <int>,
  "n_edges": <int>,
  "options": {
    "R": 6371000.0,
    "level": <int>,
    "loader": {
      "path": "builtin://icosahedral/<level>",
      "reader": "builtin_icosahedral",
      "check": "strict"
    }
  },
  "schema_version": "1.0.0"
}
```

`options.loader.path` uses the `builtin://icosahedral/<level>` scheme
rather than a filesystem path — DUO's subdivision geometry is pure math,
so no mesh files need to ship with the repo. A future `.duo` file reader
will land with the EarthSciSerialization file-format spec
(`options.loader.reader = "duo_mesh"`).

## Provenance

Each fixture was produced by the Julia reference binding via:

```julia
using EarthSciDiscretizations, JSON
using EarthSciDiscretizations: to_esm
g = EarthSciDiscretizations.grids.duo(
    loader = (
        path = "builtin://icosahedral/<level>",
        reader = "builtin_icosahedral",
        check = "strict",
    ),
    R = 6.371e6, dtype = Float64, ghosts = 0,
)
d = to_esm(g)
delete!(d, "provenance")
JSON.json(d, 2)
```

The `provenance` block is intentionally **stripped** from the committed
fixture: it is binding/platform/runtime-specific (per `GRIDS_API.md`
§6.4) and would prevent byte-identity across the four bindings. The
remainder of the payload is binding-neutral and is what the ESD walker
(see `dsc-usq`) gates on at Layer A.

JSON formatting is `JSON.print(d, 2)` (Julia `JSON.jl`'s 2-space indent,
keys sorted lexicographically) plus a trailing newline, matching the
§5.4.6 canonical wire form modulo the provenance strip described above.

## Verification

The tests in `test/test_duo_fixtures.jl` (item-tagged `:grid, :duo,
:fixtures`) load each fixture, reconstruct the grid via the Julia
binding's `duo(...)` generator, re-lower with `to_esm`, and:

- byte-diff the round-tripped declaration against the committed `.esm`
  (Layer A in `dsc-usq` walker parlance);
- check the closed-sphere Euler identity (`V - E + F = 2`);
- check topology and metric accessors (`cell_centers`, `neighbors`,
  `metric_eval`, `total_area`) at representative cells.

Other bindings (Python, Rust, TypeScript) consume these same fixtures via
their own conformance harnesses (`tests/conformance/grids/duo/`); the
binding-neutral payload guarantees the comparison is meaningful.

## Regenerating

If a deliberate change to the `duo` `.esm` lowering lands, regenerate all
fixtures from the Julia binding (the §4.3 reference binding) and commit
the diff:

```bash
julia discretizations/grids/duo/regenerate_fixtures.jl
```

The script activates a temporary env so it runs cleanly without touching
the main `Project.toml` deps.
