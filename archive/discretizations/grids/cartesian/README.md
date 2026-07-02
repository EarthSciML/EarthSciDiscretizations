# Cartesian grid fixtures

Canonical `.esm` declarative configurations for the `cartesian` grid family,
per the cross-binding API contract in `docs/GRIDS_API.md`.

Per the bead `dsc-suu` mayor correction (2026-04-20): an `.esm` grid is a
small declarative config (family + dimensions + extents + connectivity refs),
NOT a serialized geometry blob. Each binding's runtime provides accessors
(`cell_centers`, `cell_widths`, `neighbors`, `metric_eval`) that DERIVE
geometry from the declaration via pure math.

## Fixtures

| File                          | ndim | n         | extent                                | uniform | Purpose                          |
|-------------------------------|------|-----------|---------------------------------------|---------|----------------------------------|
| `uniform_1d_n16.esm`          |   1  | 16        | [0, 1]                                | yes     | Smallest canonical 1D resolution |
| `uniform_1d_n64.esm`          |   1  | 64        | [0, 1]                                | yes     | Mid 1D resolution                |
| `uniform_1d_n256.esm`         |   1  | 256       | [0, 1]                                | yes     | Largest canonical 1D resolution  |
| `uniform_2d_n64.esm`          |   2  | 64×64     | [0, 1] × [0, 1]                       | yes     | 2D unit square                   |
| `uniform_3d_n16.esm`          |   3  | 16×16×16  | [0, 1] × [0, 1] × [0, 1]              | yes     | 3D unit cube                     |
| `nonuniform_2d_stretched.esm` |   2  | 4×3       | [0, 1] × [-1, 1] (explicit edges)     | no      | Exercises the non-uniform path   |

The 1D resolutions match the bead's `N=16/64/256` ladder. The 2D / 3D / non-
uniform fixtures cover the remaining branches of `cartesian`'s constructor.

## Provenance

Each fixture was produced by the Julia binding via:

```julia
using EarthSciDiscretizations
g = EarthSciDiscretizations.grids.cartesian(; <opts shown above>)
d = to_esm(g)
delete!(d, "provenance")
JSON.json(d, 2)
```

The `provenance` block is intentionally **stripped** from the committed
fixture: it is binding/platform/runtime-specific (per `GRIDS_API.md` §6.4)
and would prevent byte-identity across the four bindings. The remainder of
the payload is binding-neutral and is what cross-binding conformance gates
on.

JSON formatting is `JSON.json(d, 2)` (Julia `JSON.jl`'s 2-space indent,
keys sorted lexicographically), matching the §5.4.6 canonical wire form
modulo the provenance strip described above.

## Verification

The tests in `test/test_cartesian_fixtures.jl` (item-tagged
`:cartesian, :fixtures`) load each fixture, reconstruct the grid via the
Julia binding's `cartesian(...)` generator, re-lower with `to_esm`, and:

- byte-diff the round-tripped declaration against the committed `.esm`
  (Layer A in `dsc-usq` walker parlance);
- spot-check topology and metric accessors (`cell_centers`, `cell_widths`,
  `cell_volume`, `neighbors`, `metric_eval`) at representative cells.

Other bindings (Python, Rust, TypeScript) consume these same fixtures via
their own conformance harnesses; the binding-neutral payload guarantees
the comparison is meaningful.

## Regenerating

If a deliberate change to the `cartesian` `.esm` lowering lands, regenerate
all fixtures from the Julia binding (the §4.3 reference binding) and commit
the diff:

```bash
julia --project=. -e '
  using EarthSciDiscretizations, JSON
  using EarthSciDiscretizations: to_esm
  for (fname, opts) in [
      ("uniform_1d_n16",  (nx=16, extent=[(0.0, 1.0)])),
      ("uniform_1d_n64",  (nx=64, extent=[(0.0, 1.0)])),
      ("uniform_1d_n256", (nx=256, extent=[(0.0, 1.0)])),
      ("uniform_2d_n64",  (nx=64, ny=64,
                           extent=[(0.0, 1.0), (0.0, 1.0)])),
      ("uniform_3d_n16",  (nx=16, ny=16, nz=16,
                           extent=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)])),
      ("nonuniform_2d_stretched", (edges=[[0.0, 0.1, 0.3, 0.7, 1.0],
                                          [-1.0, 0.0, 0.5, 1.0]],)),
  ]
      g = EarthSciDiscretizations.grids.cartesian(; opts...)
      d = to_esm(g); delete!(d, "provenance")
      open(joinpath("discretizations/grids/cartesian", fname * ".esm"), "w") do io
          JSON.print(io, d, 2); print(io, "\n")
      end
  end
'
```

All four bindings must then re-verify against the regenerated fixtures.
