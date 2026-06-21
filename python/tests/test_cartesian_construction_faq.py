"""Cross-language conformance for the cartesian **construction FAQ** (esd-dru).

Builds each fixture grid imperatively, re-derives its construction arrays through
the elementwise FAQ bridge (:func:`cartesian_construction_faq`, routed via the
ESS ``eval_coeff`` evaluator), and asserts the result is

1. byte/ULP-identical to the imperative ``CartesianGrid`` accessors (parity),
2. equal to the committed Julia-reference
   ``tests/conformance/grids/cartesian/construction/golden.json`` (the
   cross-binding contract — floats within the fixture tolerance, the integer
   neighbor maps / boundary masks exactly).

Golden arrays are stored as compact JSON strings (the dense byte form shared with
``tests/conformance/grids/duo/topology/golden.json``); neighbor ids are 0-based
with a ``-1`` off-grid sentinel and boundary masks are ``0``/``1``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_discretizations import grids
from earthsci_discretizations.grids.cartesian_faq import cartesian_construction_faq

HARNESS_DIR = (
    Path(__file__).resolve().parents[2] / "tests" / "conformance" / "grids" / "cartesian"
)
_SPEC = json.loads((HARNESS_DIR / "fixtures.json").read_text())
_GOLDEN = json.loads((HARNESS_DIR / "construction" / "golden.json").read_text())
REL_TOL = float(_GOLDEN["tolerance"]["relative"])
_BY_NAME = {f["name"]: f for f in _GOLDEN["fixtures"]}
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


def _unflatten(k: int, n: tuple[int, ...], strides: list[int]) -> tuple[int, ...]:
    return tuple((k // strides[d]) % n[d] for d in range(len(n)))


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_cartesian_construction_faq(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture["opts"])
    faq = cartesian_construction_faq(grid)
    gl = _BY_NAME[fixture_name]

    ndim = grid.ndim
    n = tuple(int(x) for x in grid.n)
    nc = 1
    for nd in n:
        nc *= nd
    strides = [1] * ndim
    for d in range(1, ndim):
        strides[d] = strides[d - 1] * n[d - 1]

    # --- counts ---
    assert ndim == int(gl["ndim"])
    assert nc == int(gl["n_cells"])

    # --- parity: FAQ == imperative accessors ---
    for d in range(ndim):
        for a, b in zip(faq.edges[d], grid.edges[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: edges[{d}] parity"
        for a, b in zip(faq.centers[d], grid.centers[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: centers[{d}] parity"
        for a, b in zip(faq.widths[d], grid.widths[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: widths[{d}] parity"
    for k in range(nc):
        idx = _unflatten(k, n, strides)
        cc = grid.cell_centers(*idx)
        cw = grid.cell_widths(*idx)
        for d in range(ndim):
            assert _close_rel(float(faq.cell_centers[d][k]), float(cc[d]))
            assert _close_rel(float(faq.cell_widths[d][k]), float(cw[d]))
        assert _close_rel(float(faq.cell_volume[k]), float(grid.cell_volume(*idx)))

    # --- byte/ULP identity against the cross-binding golden ---
    g_edges = json.loads(gl["edges"])
    g_centers = json.loads(gl["centers"])
    g_widths = json.loads(gl["widths"])
    for d in range(ndim):
        for a, b in zip(faq.edges[d], g_edges[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: edges[{d}] golden"
        for a, b in zip(faq.centers[d], g_centers[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: centers[{d}] golden"
        for a, b in zip(faq.widths[d], g_widths[d], strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}: widths[{d}] golden"

    g_vol = json.loads(gl["cell_volume"])
    g_jac = json.loads(gl["metric_jacobian"])
    for a, b in zip(faq.cell_volume, g_vol, strict=True):
        assert _close_rel(float(a), float(b)), f"{fixture_name}: cell_volume golden"
    for a, b in zip(faq.metric_jacobian, g_jac, strict=True):
        assert _close_rel(float(a), float(b)), f"{fixture_name}: metric_jacobian golden"

    g_nm = json.loads(gl["neighbor_minus"])
    g_np = json.loads(gl["neighbor_plus"])
    g_bl = json.loads(gl["boundary_lower"])
    g_bu = json.loads(gl["boundary_upper"])
    for d in range(ndim):
        assert [int(x) for x in faq.neighbor_minus[d]] == [int(x) for x in g_nm[d]]
        assert [int(x) for x in faq.neighbor_plus[d]] == [int(x) for x in g_np[d]]
        assert [int(x) for x in faq.boundary_lower[d]] == [int(x) for x in g_bl[d]]
        assert [int(x) for x in faq.boundary_upper[d]] == [int(x) for x in g_bu[d]]
