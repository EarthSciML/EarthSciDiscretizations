"""Cross-language conformance for the lat_lon grid family.

Verifies that the Python accessor runtime matches the committed reference
golden within the `docs/GRIDS_API.md` §4.2 ULP tolerance.

The Julia binding for lat_lon is tracked by bead dsc-1ts and not yet on main;
the reference golden is produced by the Python binding for this initial corpus
(see `tests/conformance/grids/latlon/regenerate_golden.py`). Python is a
self-check here — the real cross-language value comes from the Rust and
TypeScript suites loading the same golden.
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
    / "latlon"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
GOLDEN_DIR = HARNESS_DIR / "golden"
_SPEC = json.loads(FIXTURES_PATH.read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]

_METRIC_NAMES = (
    "J",
    "g_lonlon",
    "g_latlat",
    "g_lonlat",
    "ginv_lonlon",
    "ginv_latlat",
    "ginv_lonlat",
)
_DIR_KEYS = ("W", "E", "S", "N")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _build_grid(opts: dict):
    kwargs = {
        "variant": opts["variant"],
        "R": float(opts["R"]),
        "ghosts": int(opts["ghosts"]),
        "dtype": opts["dtype"],
        "pole_policy": opts["pole_policy"],
    }
    if opts["variant"] == "regular":
        kwargs["nlon"] = int(opts["nlon"])
        kwargs["nlat"] = int(opts["nlat"])
    else:
        kwargs["nlon_per_row"] = [int(x) for x in opts["nlon_per_row"]]
        if "lat_edges" in opts:
            kwargs["lat_edges"] = [float(x) for x in opts["lat_edges"]]
    return grids.lat_lon(**kwargs)


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture["opts"])
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    assert grid.n_cells == golden["n_cells"]
    assert len(fixture["query_points"]) == len(golden["cell_centers"])

    for k, (j, i) in enumerate(fixture["query_points"]):
        lon, lat = grid.cell_centers(j=j, i=i)
        gc = golden["cell_centers"][k]
        assert _close_rel(lon, gc["lon"]), f"lon at qp[{k}]=({j},{i})"
        assert _close_rel(lat, gc["lat"]), f"lat at qp[{k}]=({j},{i})"

        nbrs = grid.neighbors(j, i)
        for dkey in _DIR_KEYS:
            expected = golden[f"neighbors_{dkey}"][k]
            got = nbrs[dkey]
            if expected is None:
                assert got is None, (
                    f"neighbor {dkey} at qp[{k}]=({j},{i}): expected pole, got {got}"
                )
            else:
                assert got is not None, (
                    f"neighbor {dkey} at qp[{k}]=({j},{i}): expected {expected}, got None"
                )
                got_list = [int(got[0]), int(got[1])]
                exp_list = [int(expected[0]), int(expected[1])]
                assert got_list == exp_list, (
                    f"neighbor {dkey} at qp[{k}]=({j},{i}): got {got_list}, expected {exp_list}"
                )

        for name in _METRIC_NAMES:
            got = grid.metric_eval(name, j, i)
            expected = float(golden[f"metric_{name}"][k])
            assert _close_rel(got, expected), (
                f"metric {name} at qp[{k}]=({j},{i}): got {got}, expected {expected}"
            )

        got_area = grid.metric_eval("area", j, i)
        expected_area = float(golden["area"][k])
        assert _close_rel(got_area, expected_area), (
            f"area at qp[{k}]=({j},{i}): got {got_area}, expected {expected_area}"
        )
