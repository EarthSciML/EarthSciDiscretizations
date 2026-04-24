# Canonical cubed_sphere grid fixtures

Declarative `.esm` grid instances for the `cubed_sphere` family at standard
resolutions. Per the 2026-04-20 mayor correction, each `.esm` is a small
declarative config (family + dimensions + extents + generator reference),
**not** a serialized geometry blob. Accessor outputs are derived at runtime
by the binding from the declaration (see `docs/GRIDS_API.md` §6).

## Fixtures

| File        | Nc  | R (m)     | Notes                                          |
|-------------|-----|-----------|------------------------------------------------|
| `c24.esm`   | 24  | 6.371e6   | Coarse (~370 km). Fast smoke tests.            |
| `c48.esm`   | 48  | 6.371e6   | Typical research resolution (~180 km).         |
| `c96.esm`   | 96  | 6.371e6   | Production climate resolution (~90 km).        |

All fixtures use `dtype: float64`, `ghosts: 0`, and the builtin `gnomonic_c6`
generator (equiangular gnomonic projection on a 6-panel cube).

## Schema shape

Each fixture conforms to the §6 declarative grid schema as lowered by the
binding `to_esm` contract (see `python/.../cubed_sphere.py::to_esm` and
`rust/.../cubed_sphere.rs::Grid::to_esm` for reference implementations):

```json
{
  "family": "cubed_sphere",
  "version": "1.0.0",
  "dtype": "float64",
  "topology": "block_structured",
  "generator": "gnomonic_c6",
  "params": { "Nc": <int>, "R": <float>, "ghosts": <int> }
}
```

Inline runtime `to_esm` also emits a `provenance` block (§6.4) that records
the emitting binding, version, platform, and math-lib fingerprint. The
committed fixtures **omit** `provenance` because it is per-build state;
cross-binding byte-identity is promised on the declarative portion alone
(per §6.3, provenance is normalized by `canonicalize()`).

## Provenance

- **Reference binding**: Julia (per `docs/GRIDS_API.md` §4.3 tie-break rule).
- **Generated via**: the canonical declarative shape emitted by
  `EarthSciDiscretizations.grids.cubed_sphere(Nc = N, R = 6.371e6).to_esm()`
  with `provenance` stripped.
- **Verified cross-binding**: the same `(family, generator, params)` tuple
  produces bit-identical accessor outputs under the cubed-sphere conformance
  harness at `tests/conformance/grids/cubed_sphere/` (landed in `dsc-fgq`),
  within the tolerance published in `fixtures.json`.

## Consumed by

- `test/test_cubed_sphere_fixtures.jl` — Layer A byte-diff checks plus
  topology/metric sanity asserts against each fixture (consumed by the ESD
  CI walker `dsc-usq`).
- Cross-binding conformance harness at
  `tests/conformance/grids/cubed_sphere/fixtures.json` — uses the same
  `(Nc, R)` tuples at the accessor level.

## Regenerating

If the declarative schema changes (§6 minor or major bump), regenerate:

```bash
julia --project=. -e '
  using EarthSciDiscretizations, JSON
  for (name, Nc) in (("c24", 24), ("c48", 48), ("c96", 96))
      doc = Dict(
          "family"    => "cubed_sphere",
          "version"   => "1.0.0",
          "dtype"     => "float64",
          "topology"  => "block_structured",
          "generator" => "gnomonic_c6",
          "params"    => Dict("Nc" => Nc, "R" => 6.371e6, "ghosts" => 0),
      )
      open("discretizations/grids/cubed_sphere/$name.esm", "w") do io
          JSON.print(io, doc, 2)
          println(io)
      end
  end
'
```

Bump `version` in lockstep with the schema bump and re-verify all four
bindings against the new corpus.
