# earthsci-toolkit (Python binding for EarthSciDiscretizations)

Python implementation of the EarthSciDiscretizations grid accessors and
discretization runtime. Conforms to the cross-binding contract in
[`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Install

```bash
pip install -e .[dev]
# Optional ecosystem integration:
pip install -e .[dev,xarray]
```

## Layout

- `src/earthsci_toolkit/` — package source
- `src/earthsci_toolkit/grids/` — per-family grid accessor modules (one file per family)
- `tests/` — pytest suite

## Test

```bash
pytest
```

See `docs/GRIDS_API.md` (repo root) for the normative API contract this
package implements.
