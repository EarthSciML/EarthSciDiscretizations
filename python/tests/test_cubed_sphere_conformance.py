"""Cross-language conformance for the cubed_sphere grid family.

Verifies that the Python accessor runtime produces accessor outputs that
match the committed Julia-reference golden within the documented
`docs/GRIDS_API.md` §4.2 ULP tolerance. Corpus lives at
``tests/conformance/grids/cubed_sphere/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_toolkit import grids

HARNESS_DIR = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "conformance"
    / "grids"
    / "cubed_sphere"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
GOLDEN_DIR = HARNESS_DIR / "golden"
_SPEC = json.loads(FIXTURES_PATH.read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    opts = fixture["opts"]
    grid = grids.cubed_sphere(
        Nc=opts["Nc"], R=opts["R"], ghosts=opts["ghosts"], dtype=opts["dtype"]
    )
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    assert grid.n_cells == golden["n_cells"]
    assert len(fixture["query_points"]) == len(golden["cell_centers"])

    metric_names = (
        "J",
        "g_xixi",
        "g_etaeta",
        "g_xieta",
        "ginv_xixi",
        "ginv_etaeta",
        "ginv_xieta",
    )

    for k, (p, i, j) in enumerate(fixture["query_points"]):
        lon, lat = grid.cell_centers(p=p, i=i, j=j)
        gc = golden["cell_centers"][k]
        assert _close_rel(lon, gc["lon"]), f"lon mismatch at qp[{k}]=({p},{i},{j})"
        assert _close_rel(lat, gc["lat"]), f"lat mismatch at qp[{k}]=({p},{i},{j})"

        nbrs = grid.neighbors(p, i, j)
        for dkey in ("W", "E", "S", "N"):
            got = list(nbrs[dkey])
            expected = [int(x) for x in golden[f"neighbors_{dkey}"][k]]
            assert got == expected, f"neighbor {dkey} mismatch at qp[{k}]=({p},{i},{j})"

        for name in metric_names:
            got = grid.metric_eval(name, p, i, j)
            expected = float(golden[f"metric_{name}"][k])
            assert _close_rel(got, expected), (
                f"metric {name} mismatch at qp[{k}]=({p},{i},{j}): got {got}, expected {expected}"
            )

        got_area = grid.metric_eval("area", p, i, j)
        expected_area = float(golden["area"][k])
        assert _close_rel(got_area, expected_area), (
            f"area mismatch at qp[{k}]=({p},{i},{j})"
        )
