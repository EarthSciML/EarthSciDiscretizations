# Canonical arakawa grid instances

Canonical `.esm` declarative configurations for the `arakawa` grid family at
standard resolutions. Committed per bead `dsc-u0l`:
`[Grids/arakawa/fixtures] Canonical arakawa grid instances in ESD`.

## Scope

Per the 2026-04-20 mayor correction (bead `dsc-suu`): a `.esm` grid entry is a
**small declarative config** (family + base-grid ref + stagger + dimensions +
extents), **not** a serialized geometry blob. The runtime in each binding
*derives* geometry on demand via accessors (`cell_centers`, `neighbors`,
`metric_eval`). Cross-binding conformance therefore compares accessor outputs
at pinned query points (see `tests/conformance/grids/cubed_sphere/` for the
per-family pattern), not serialized bytes.

These fixtures are the canonical configuration inputs a future arakawa
conformance harness would feed to each binding's generator.

## Contents

| File                          | Base   | Resolution           | Stagger |
|-------------------------------|--------|----------------------|---------|
| `cartesian_n16.esm.json`      | cart.  | nx=ny=16 on [0,1]²   | C       |
| `cartesian_n64.esm.json`      | cart.  | nx=ny=64 on [0,1]²   | C       |
| `cartesian_n256.esm.json`     | cart.  | nx=ny=256 on [0,1]²  | C       |
| `lat_lon_1deg.esm.json`       | latlon | 360 × 180 (1°)       | C       |
| `lat_lon_0p25deg.esm.json`    | latlon | 1440 × 720 (0.25°)   | C       |
| `lat_lon_0p1deg.esm.json`     | latlon | 3600 × 1800 (0.1°)   | C       |

All instances use `stagger = "C"` (cell-centered scalars; u on E/W faces, v on
N/S faces) — the dominant staggering in atmosphere/ocean finite-volume
dynamical cores. Other staggers (A/B/D/E) are verified by the inline
`test/test_arakawa.jl` unit tests, not the canonical fixture set.

`fixtures.json` is the manifest consumed by the inline conformance tests.

## Canonical form

Every `.esm.json` contains the declarative fields:

- `family: "arakawa"`
- `stagger: "A" | "B" | "C" | "D" | "E"`
- `rotated: bool` (true only when `stagger == "E"`)
- `topology: "block_structured"`
- `dtype: "float64" | "float32"`
- `ghosts: int ≥ 0`
- `n_cells: int` (= `base.nx * base.ny` or `base.nlon * base.nlat`)
- `base: {…}` — declarative ref to the base grid (see below)

Keys are sorted alphabetically within each object (canonical JSON form; the
canonicalization spec lives in EarthSciSerialization §5.4.6).

### Cartesian base ref

Matches `../cartesian.schema.json` required fields:

```json
{ "family": "cartesian", "nx": <int>, "ny": <int>, "extent": [[xlo, ylo], [xhi, yhi]] }
```

### Lat-lon base ref

Matches `../lat_lon.schema.json` required/common fields:

```json
{
  "family": "lat_lon", "variant": "regular",
  "nlon": <int>, "nlat": <int>, "R": <float>,
  "lon_start": <float>, "pole_policy": "none"
}
```

## What is not in the canonical form

- `provenance` (binding/platform/sha metadata) — **excluded**. Each binding's
  runtime adds a `provenance` block when *emitting*, but it is not part of
  the canonical declarative form used for comparison. See docs/GRIDS_API.md
  §6.4 for the provenance contract.
- `version` — omitted pending cross-binding alignment (Julia does not emit
  `version`; Python/Rust/TS do). When the Julia runtime adopts `version`,
  these fixtures regenerate to carry it.
- Inline geometry arrays (`cells`, `edges`, `vertices`) — forbidden per the
  2026-04-20 mayor correction.

## Lat-lon fixtures: declaration-only

The arakawa accessor runtime is parameterized over an `ArakawaBaseGrid`
abstraction. As of this bead (2026-04-24), all four bindings implement only
`CartesianBase`; standalone `lat_lon` generators exist (Phase 3), but an
`ArakawaLatLonBase` wrapper that lets `arakawa(base=..., stagger=...)` take a
lat-lon base has not landed. The three `lat_lon_*.esm.json` fixtures are
therefore **declaration-only** and document the intended canonical shape for
the follow-up composition bead. The inline test validates their structure
(schema-level) but does not round-trip them through a generator.

## Regenerating

The cartesian fixtures are produced by the Julia reference binding (per
`docs/GRIDS_API.md` §4.3). To regenerate after a deliberate algorithm change:

```julia
using EarthSciDiscretizations
using JSON
for (name, nx, ny) in (("cartesian_n16", 16, 16),
                       ("cartesian_n64", 64, 64),
                       ("cartesian_n256", 256, 256))
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = nx, ny = ny)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)
    open("$(name).esm.json", "w") do io
        JSON.print(io, to_esm(g), 2)
    end
end
```

Commit the regenerated files; the inline test (`test/test_arakawa_fixtures.jl`)
auto-verifies them on the next CI run.
