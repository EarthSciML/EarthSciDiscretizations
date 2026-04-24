"""Cross-language conformance for the cartesian grid family.

Verifies that the Python accessor runtime produces accessor outputs that
match the committed Julia-reference golden within the documented
``docs/GRIDS_API.md`` §4.2 tolerance. Corpus lives at
``tests/conformance/grids/cartesian/``.
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
    / "cartesian"
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


def _build_grid(opts: dict):
    kw: dict = {"dtype": opts.get("dtype", "float64"), "ghosts": int(opts.get("ghosts", 0))}
    if "edges" in opts:
        kw["edges"] = [list(map(float, e)) for e in opts["edges"]]
    else:
        if "nx" in opts:
            kw["nx"] = int(opts["nx"])
        if "ny" in opts:
            kw["ny"] = int(opts["ny"])
        if "nz" in opts:
            kw["nz"] = int(opts["nz"])
        if "extent" in opts:
            kw["extent"] = [(float(lo), float(hi)) for lo, hi in opts["extent"]]
    return grids.cartesian(**kw)


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture["opts"])
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    ndim = grid.ndim
    assert ndim == int(golden["ndim"])
    assert grid.n_cells == int(golden["n_cells"])
    assert len(fixture["query_points"]) == len(golden["cell_centers"])

    width_names = ("dx", "dy", "dz")
    face_names = ("face_area_x", "face_area_y", "face_area_z")

    for k, qp in enumerate(fixture["query_points"]):
        idx = tuple(int(x) for x in qp)

        c = grid.cell_centers(*idx)
        gc = golden["cell_centers"][k]
        for d in range(ndim):
            assert _close_rel(float(c[d]), float(gc[d])), (
                f"{fixture_name}: cell_centers[{d}] at qp[{k}]={idx}"
            )

        w = grid.cell_widths(*idx)
        gw = golden["cell_widths"][k]
        for d in range(ndim):
            assert _close_rel(float(w[d]), float(gw[d])), (
                f"{fixture_name}: cell_widths[{d}] at qp[{k}]={idx}"
            )

        v = grid.cell_volume(*idx)
        assert _close_rel(float(v), float(golden["cell_volume"][k])), (
            f"{fixture_name}: cell_volume at qp[{k}]={idx}"
        )

        # Neighbors: Python returns a dict keyed by (axis, side) with 0-based
        # axis indices. Serialize to the same sorted list-of-entries shape
        # the golden emits (axis-major, low side before high side).
        nbr_dict = grid.neighbors(*idx)
        nbr_list = []
        for d in range(ndim):
            for s in (-1, +1):
                if (d, s) in nbr_dict:
                    nbr_list.append({
                        "axis": d,
                        "side": s,
                        "index": list(int(x) for x in nbr_dict[(d, s)]),
                    })
        gnbrs = golden["neighbors"][k]
        assert len(nbr_list) == len(gnbrs), (
            f"{fixture_name}: neighbor count at qp[{k}]={idx}"
        )
        for got, exp in zip(nbr_list, gnbrs):
            assert got["axis"] == int(exp["axis"])
            assert got["side"] == int(exp["side"])
            assert got["index"] == [int(x) for x in exp["index"]]

        # Scalar metrics
        assert _close_rel(
            float(grid.metric_eval("volume", *idx)),
            float(golden["metric_volume"][k]),
        )
        assert _close_rel(
            float(grid.metric_eval("jacobian", *idx)),
            float(golden["metric_jacobian"][k]),
        )
        for d in range(ndim):
            assert _close_rel(
                float(grid.metric_eval(width_names[d], *idx)),
                float(golden[f"metric_{width_names[d]}"][k]),
            ), f"{fixture_name}: metric {width_names[d]} at qp[{k}]={idx}"
            assert _close_rel(
                float(grid.metric_eval(face_names[d], *idx)),
                float(golden[f"metric_{face_names[d]}"][k]),
            ), f"{fixture_name}: metric {face_names[d]} at qp[{k}]={idx}"

        g = grid.metric_eval("g", *idx)
        gg = golden["metric_g"][k]
        for i in range(ndim):
            for j in range(ndim):
                assert _close_rel(float(g[i][j]), float(gg[i][j])), (
                    f"{fixture_name}: metric g[{i}][{j}] at qp[{k}]={idx}"
                )
