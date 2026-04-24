"""Cross-language conformance for the mpas grid family.

Verifies that the Python accessor runtime produces accessor outputs that
match the committed reference golden within the documented
``docs/GRIDS_API.md`` §4.2 ULP tolerance. Corpus lives at
``tests/conformance/grids/mpas/``.

The mpas reference binding is interim Python until dsc-7j0 (Julia mpas)
lands — see ``tests/conformance/grids/mpas/README.md`` for the rationale.
Regardless of which binding generated the golden, every binding's
accessors must match it within the declared tolerance at the pinned
query points.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from earthsci_toolkit import grids
from earthsci_toolkit.grids import mpas_mesh_data

HARNESS_DIR = (
    Path(__file__).resolve().parents[2] / "tests" / "conformance" / "grids" / "mpas"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
GOLDEN_DIR = HARNESS_DIR / "golden"
_SPEC = json.loads(FIXTURES_PATH.read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]

_CELL_METRICS = ("lon", "lat", "area", "x", "y", "z", "n_edges_on_cell")
_EDGE_METRICS = ("lon_edge", "lat_edge", "dc_edge", "dv_edge")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _build_grid(fixture: dict):
    m = fixture["mesh"]
    opts = fixture["opts"]
    mesh = mpas_mesh_data(
        lon_cell=np.asarray(m["lon_cell"], dtype=np.float64),
        lat_cell=np.asarray(m["lat_cell"], dtype=np.float64),
        area_cell=np.asarray(m["area_cell"], dtype=np.float64),
        n_edges_on_cell=np.asarray(m["n_edges_on_cell"], dtype=np.int32),
        cells_on_cell=np.asarray(m["cells_on_cell"], dtype=np.int32),
        edges_on_cell=np.asarray(m["edges_on_cell"], dtype=np.int32),
        lon_edge=np.asarray(m["lon_edge"], dtype=np.float64),
        lat_edge=np.asarray(m["lat_edge"], dtype=np.float64),
        cells_on_edge=np.asarray(m["cells_on_edge"], dtype=np.int32),
        dc_edge=np.asarray(m["dc_edge"], dtype=np.float64),
        dv_edge=np.asarray(m["dv_edge"], dtype=np.float64),
        max_edges=int(m["max_edges"]),
        x_cell=np.asarray(m["x_cell"], dtype=np.float64),
        y_cell=np.asarray(m["y_cell"], dtype=np.float64),
        z_cell=np.asarray(m["z_cell"], dtype=np.float64),
        n_vertices=int(m["n_vertices"]),
        R=float(opts["R"]),
    )
    return grids.mpas(
        mesh=mesh,
        R=float(opts["R"]),
        dtype=str(opts["dtype"]),
        ghosts=int(opts["ghosts"]),
    )


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture)
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    assert grid.n_cells == golden["n_cells"]
    assert grid.n_edges == golden["n_edges"]
    assert grid.max_edges == golden["max_edges"]

    for k, c in enumerate(fixture["query_cells"]):
        lon, lat = grid.cell_centers(c)
        gc = golden["cell_centers"][k]
        assert _close_rel(lon, gc["lon"]), f"lon mismatch at cell {c}"
        assert _close_rel(lat, gc["lat"]), f"lat mismatch at cell {c}"

        x, y, z = grid.cell_center_cart(c)
        gx = golden["cell_centers_cart"][k]
        assert _close_rel(x, gx["x"]), f"cart x mismatch at cell {c}"
        assert _close_rel(y, gx["y"]), f"cart y mismatch at cell {c}"
        assert _close_rel(z, gx["z"]), f"cart z mismatch at cell {c}"

        got_nb = [int(n) for n in grid.neighbors(c)]
        exp_nb = [int(n) for n in golden["neighbors"][k]]
        assert got_nb == exp_nb, f"neighbors mismatch at cell {c}"

        assert _close_rel(grid.cell_area(c), float(golden["cell_area"][k])), (
            f"cell_area mismatch at cell {c}"
        )

        for name in _CELL_METRICS:
            got = grid.metric_eval(name, c)
            expected = float(golden["cell_metrics"][name][k])
            assert _close_rel(got, expected), f"cell metric {name} mismatch at cell {c}"

    for k, e in enumerate(fixture["query_edges"]):
        assert _close_rel(grid.edge_length(e), float(golden["edge_length"][k])), (
            f"edge_length mismatch at edge {e}"
        )
        assert _close_rel(grid.cell_distance(e), float(golden["cell_distance"][k])), (
            f"cell_distance mismatch at edge {e}"
        )
        for name in _EDGE_METRICS:
            got = grid.metric_eval(name, e)
            expected = float(golden["edge_metrics"][name][k])
            assert _close_rel(got, expected), f"edge metric {name} mismatch at edge {e}"
