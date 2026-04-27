"""Cross-language conformance for the DUO grid family.

Verifies that the Python accessor runtime produces accessor outputs that
match the committed Julia-reference golden within the documented
``docs/GRIDS_API.md`` §4.2 ULP tolerance. Corpus lives at
``tests/conformance/grids/duo/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_discretizations import grids

HARNESS_DIR = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "conformance"
    / "grids"
    / "duo"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
GOLDEN_DIR = HARNESS_DIR / "golden"
_SPEC = json.loads(FIXTURES_PATH.read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]

_METRIC_NAMES = ("area", "lon", "lat", "x", "y", "z")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_duo_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    opts = fixture["opts"]
    level = int(opts["level"])
    grid = grids.duo(
        loader={"path": f"builtin://icosahedral/{level}"},
        R=float(opts["R"]),
        ghosts=int(opts["ghosts"]),
        dtype=opts["dtype"],
    )
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    assert grid.n_cells == golden["n_cells"]
    assert grid.n_vertices == golden["n_vertices"]
    assert grid.n_edges == golden["n_edges"]
    assert len(fixture["query_points"]) == len(golden["cell_centers"])

    for k, c in enumerate(fixture["query_points"]):
        c = int(c)
        lon, lat = grid.cell_centers(c)
        gc = golden["cell_centers"][k]
        assert _close_rel(lon, gc["lon"]), f"lon mismatch at qp[{k}]=c{c}"
        assert _close_rel(lat, gc["lat"]), f"lat mismatch at qp[{k}]=c{c}"

        got_nb = list(grid.neighbors(c))
        expected_nb = [int(x) for x in golden["neighbors"][k]]
        assert got_nb == expected_nb, f"neighbors mismatch at qp[{k}]=c{c}"

        for m in _METRIC_NAMES:
            got = grid.metric_eval(m, c)
            expected = float(golden[f"metric_{m}"][k])
            assert _close_rel(got, expected), (
                f"metric {m} mismatch at qp[{k}]=c{c}: got {got}, expected {expected}"
            )
