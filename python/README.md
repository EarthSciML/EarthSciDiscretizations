# earthsci-discretizations (Python binding for EarthSciDiscretizations)

Python implementation of the EarthSciDiscretizations grid accessors and
discretization runtime. Conforms to the cross-binding contract in
[`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Install

This package depends on the unpublished `earthsci-toolkit` from
[EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization);
pip resolves it via a `git+` URL with `subdirectory=packages/earthsci_toolkit`.

```bash
pip install -e .[dev]
# The `xarray` extra (xarray + pyproj) is declared for the planned
# ecosystem grid views (docs/GRIDS_API.md §5.2); no integration code
# consumes it yet:
pip install -e .[dev,xarray]
```

## Layout

- `src/earthsci_discretizations/` — package source
- `src/earthsci_discretizations/grids/` — per-family grid accessor modules (one file per family)
- `src/earthsci_discretizations/rules/` — discretization-rule runtime (AST
  `eval_coeff` adapter over `earthsci_toolkit.evaluate`, plus the
  `load_rule` / `Rule` / `StencilEntry` rule loader). Stencil application,
  ghost-cell synthesis, and parabola reconstruction are owned by the
  canonical `earthsci_toolkit` pipeline, not duplicated here. Mirrors
  `src/rule_eval.jl` in the Julia binding.
- `tests/` — pytest suite

## Test

```bash
pytest
```

See `docs/GRIDS_API.md` (repo root) for the normative API contract this
package implements.
