"""Cross-language conformance for the vertical **construction FAQ** (esd-dru).

Builds each fixture grid imperatively, re-derives its construction arrays through
the elementwise FAQ bridge (:func:`vertical_construction_faq`, routed via the ESS
``eval_coeff`` evaluator), and asserts the result is

1. byte-identical to the imperative ``VerticalGrid`` accessors (parity),
2. equal to the committed Julia-reference
   ``tests/conformance/grids/vertical/construction/golden.json``.

The vertical golden carries ``tolerance.relative = 0.0`` (no transcendentals —
strict byte equality). Arrays are stored as compact JSON strings; neighbor ids
are 0-based with a ``-1`` off-column sentinel and boundary masks are 0/1.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_discretizations import grids
from earthsci_discretizations.grids.vertical_faq import vertical_construction_faq

HARNESS_DIR = (
    Path(__file__).resolve().parents[2] / "tests" / "conformance" / "grids" / "vertical"
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
    kw: dict = {"coordinate": opts["coordinate"]}
    for key in ("nz", "ghosts"):
        if key in opts:
            kw[key] = int(opts[key])
    for key in ("levels", "ak", "bk"):
        if key in opts:
            kw[key] = [float(x) for x in opts[key]]
    if "p0" in opts:
        kw["p0"] = float(opts["p0"])
    return grids.vertical(**kw)


def _assert_floats(got, golden_str: str, what: str) -> None:
    g = json.loads(golden_str)
    assert len(got) == len(g), f"{what}: length"
    for a, b in zip(got, g, strict=True):
        assert _close_rel(float(a), float(b)), f"{what}: {a} vs {b}"


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_vertical_construction_faq(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture["opts"])
    faq = vertical_construction_faq(grid)
    gl = _BY_NAME[fixture_name]
    n = len(faq.widths)

    # --- counts ---
    assert n == int(gl["n_cells"])
    assert len(faq.levels) == int(gl["n_vertices"])

    # --- parity: FAQ == imperative accessors ---
    for a, b in zip(faq.levels, grid.levels, strict=True):
        assert _close_rel(float(a), float(b)), f"{fixture_name}: levels parity"
    for a, b in zip(faq.centers, grid.centers, strict=True):
        assert _close_rel(float(a), float(b)), f"{fixture_name}: centers parity"
    for a, b in zip(faq.widths, grid.widths, strict=True):
        assert _close_rel(float(a), float(b)), f"{fixture_name}: widths parity"
    for k in range(n):
        assert _close_rel(float(faq.metric_dz[k]), float(grid.metric_eval("dz", k)))
        assert _close_rel(float(faq.metric_z[k]), float(grid.metric_eval("z", k)))

    # --- byte identity against the cross-binding golden ---
    _assert_floats(faq.levels, gl["levels"], "levels")
    _assert_floats(faq.centers, gl["centers"], "centers")
    _assert_floats(faq.widths, gl["widths"], "widths")
    _assert_floats(faq.cell_volume, gl["cell_volume"], "cell_volume")
    _assert_floats(faq.metric_dz, gl["metric_dz"], "metric_dz")
    _assert_floats(faq.metric_z, gl["metric_z"], "metric_z")

    # Optional metrics — present in the golden iff the coordinate defines them.
    if "metric_sigma" in gl:
        assert faq.metric_sigma is not None
        _assert_floats(faq.metric_sigma, gl["metric_sigma"], "metric_sigma")
    if "metric_pressure" in gl:
        assert faq.metric_pressure is not None
        _assert_floats(faq.metric_pressure, gl["metric_pressure"], "metric_pressure")
    if "metric_ak" in gl:
        assert faq.metric_ak is not None
        _assert_floats(faq.metric_ak, gl["metric_ak"], "metric_ak")
    if "metric_bk" in gl:
        assert faq.metric_bk is not None
        _assert_floats(faq.metric_bk, gl["metric_bk"], "metric_bk")
    if "ak" in gl:
        _assert_floats(faq.ak, gl["ak"], "ak")
    if "bk" in gl:
        _assert_floats(faq.bk, gl["bk"], "bk")
    assert _close_rel(float(faq.p0), float(gl["p0"]))

    # --- structural integer arrays (exact) ---
    assert [int(x) for x in faq.neighbor_minus] == list(json.loads(gl["neighbor_minus"]))
    assert [int(x) for x in faq.neighbor_plus] == list(json.loads(gl["neighbor_plus"]))
    assert [int(x) for x in faq.boundary_lower] == list(json.loads(gl["boundary_lower"]))
    assert [int(x) for x in faq.boundary_upper] == list(json.loads(gl["boundary_upper"]))
