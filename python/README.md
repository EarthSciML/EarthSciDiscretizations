# earthsci-discretizations (Python binding for EarthSciDiscretizations)

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

- `src/earthsci_discretizations/` — package source
- `src/earthsci_discretizations/grids/` — per-family grid accessor modules (one file per family)
- `src/earthsci_discretizations/rules/` — discretization-rule runtime (AST
  `eval_coeff`, rule loader, structured stencil application). Mirrors
  `src/rule_eval.jl` in the Julia binding.
- `tests/` — pytest suite

## Test

```bash
pytest
```

See `docs/GRIDS_API.md` (repo root) for the normative API contract this
package implements.
