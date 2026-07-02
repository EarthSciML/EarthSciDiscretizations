# Canonical lat-lon grid fixtures

Declarative `.esm` configurations for the lat-lon grid family at three
standard resolutions. Per the 2026-04-20 mayor correction (bead `dsc-brl`
NOTES), each `.esm` is a **declarative config** — family, variant,
dimensions, generator ref — not a serialized geometry blob. Geometry is
derived on demand by each binding's accessors (`cell_centers`, `neighbors`,
`metric_eval`, `cell_area`) per `docs/GRIDS_API.md` §2–§3.

Schema: [`../lat_lon.schema.json`](../lat_lon.schema.json).
API contract: [`../../../docs/GRIDS_API.md`](../../../docs/GRIDS_API.md).

## Fixtures

| File                  | nlon | nlat | Nominal resolution | Cells     |
|-----------------------|------|------|--------------------|-----------|
| `lat_lon_1deg.esm`    | 360  | 180  | 1°                 | 64 800    |
| `lat_lon_0p25deg.esm` | 1440 | 720  | 0.25°              | 1 036 800 |
| `lat_lon_0p1deg.esm`  | 3600 | 1800 | 0.1°               | 6 480 000 |

All use the `regular` variant, equal-angle latitude edges (schema default),
Earth radius `R = 6.371e6 m`, `ghosts = 0`, `pole_policy = "none"`, and the
`[-π, π)` longitude convention (`lon_start = -π`).

## Provenance

Generated from the Python binding at `python/src/earthsci_toolkit/grids/lat_lon.py`
(bead `dsc-4rj`) as the cross-binding reference for the lat-lon family:

```python
from earthsci_toolkit.grids.lat_lon import lat_lon
doc = lat_lon(nlon=360, nlat=180).to_esm()
# Equivalent calls from the Rust, TypeScript, and (pending) Julia bindings
# produce the same declarative payload under §5.4.6 canonicalization.
```

The committed fixtures omit the per-binding `provenance` block (binding
name, platform triple, math-lib fingerprint) described in
`docs/GRIDS_API.md` §6.4. That block is attached by the generator at
emit-time and is not part of the canonical declarative payload: any
binding's `to_esm` output, with `provenance` stripped, must match these
bytes under §5.4.6 canonicalization.

## Consumers

- Cross-language conformance harness (bead `dsc-zej`) —
  `tests/conformance/grids/latlon/`.
- Inline schema / topology tests —
  [`../../../test/test_latlon_fixtures.jl`](../../../test/test_latlon_fixtures.jl).
- Downstream discretization rules whose `grid_family = "lat_lon"`.

## Regenerating

```bash
python3 -c "
import sys, json
sys.path.insert(0, 'python/src')
from earthsci_toolkit.grids.lat_lon import lat_lon
for nlon, nlat, label in [(360, 180, '1deg'), (1440, 720, '0p25deg'), (3600, 1800, '0p1deg')]:
    d = lat_lon(nlon=nlon, nlat=nlat).to_esm()
    d.pop('provenance', None)
    with open(f'discretizations/grids/latlon/lat_lon_{label}.esm', 'w') as f:
        json.dump(d, f, indent=2, sort_keys=True)
        f.write('\n')
"
```

Regenerate only after a deliberate schema or generator change; any update
must keep all four bindings byte-identical under canonicalize().
