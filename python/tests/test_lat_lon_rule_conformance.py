"""Cross-binding rule-eval conformance test for the centered_2nd_uniform_latlon rule.

Verifies that ESD Python's ``eval_coeff`` (→ ``earthsci_toolkit.evaluate``) produces
coefficient values matching the Julia reference golden in
``tests/conformance/rules/centered_2nd_uniform_latlon/`` within the tolerance
specified in ``fixtures.json``.

Two levels of coverage:

1. ``coefficient_cases`` — direct bindings-to-coefficient verification.
   Each case provides explicit scalar bindings; ``eval_coeff`` is called once
   per stencil entry and compared against the Julia reference value.

2. ``stencil_cases`` — whole-stencil application on a real lat-lon grid.
   Per-cell bindings are derived from the grid accessor; the stencil sum
   applied to ``f(lon, lat) = sin(lon)*cos(lat)`` is compared against the
   Julia reference result.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from earthsci_discretizations import grids
from earthsci_discretizations.rules import eval_coeff

HARNESS_DIR = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "conformance"
    / "rules"
    / "centered_2nd_uniform_latlon"
)
REPO_ROOT = Path(__file__).resolve().parents[2]

_SPEC = json.loads((HARNESS_DIR / "fixtures.json").read_text())
_GOLDEN = json.loads((HARNESS_DIR / "golden" / "coefficients.json").read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _load_rule() -> dict:
    rule_path = REPO_ROOT / _SPEC["rule_path"]
    data = json.loads(rule_path.read_text())
    return data["discretizations"][_SPEC["rule"]]


def _selector_key(selector: dict) -> str:
    """Map a selector dict to its golden coefficient key (e.g. ``'lon_m1'``)."""
    axis = selector["axis"]
    offset = int(selector["offset"])
    suffix = "m1" if offset == -1 else "p1"
    return f"{axis}_{suffix}"


_COEFF_CASES = {c["name"]: c for c in _SPEC["coefficient_cases"]}
_GOLDEN_COEFF = {c["name"]: c for c in _GOLDEN["coefficient_cases"]}
_COEFF_CASE_NAMES = [c["name"] for c in _SPEC["coefficient_cases"]]


@pytest.mark.parametrize("case_name", _COEFF_CASE_NAMES)
def test_coefficient_case(case_name: str) -> None:
    """eval_coeff matches the Julia-reference coefficient for each binding set."""
    rule = _load_rule()
    bindings = {k: float(v) for k, v in _COEFF_CASES[case_name]["bindings"].items()}
    golden_coeffs = _GOLDEN_COEFF[case_name]["coefficients"]

    for entry in rule["stencil"]:
        key = _selector_key(entry["selector"])
        expected = float(golden_coeffs[key])
        got = eval_coeff(entry["coeff"], bindings)
        assert _close_rel(got, expected), (
            f"{case_name}: coeff '{key}' got {got}, expected {expected}"
        )


def _bindings_for(grid: grids.LatLonGrid, j: int, i: int) -> dict[str, float]:
    """Per-cell bindings required by centered_2nd_uniform_latlon."""
    _, lat = grid.cell_center(j, i)
    nlon = grid.nlon(j)
    dlon = 2.0 * math.pi / nlon
    dlat = float(grid.lat_edges[j + 1]) - float(grid.lat_edges[j])
    return {
        "R": float(grid.R),
        "cos_lat": math.cos(lat),
        "dlon": dlon,
        "dlat": dlat,
    }


_STENCIL_CASES = {c["name"]: c for c in _SPEC["stencil_cases"]}
_GOLDEN_STENCIL = {c["name"]: c for c in _GOLDEN["stencil_cases"]}
_STENCIL_CASE_NAMES = [c["name"] for c in _SPEC["stencil_cases"]]

_DIR_FOR_SELECTOR: dict[tuple[str, int], str] = {
    ("lon", -1): "W",
    ("lon", 1): "E",
    ("lat", -1): "S",
    ("lat", 1): "N",
}


@pytest.mark.parametrize("case_name", _STENCIL_CASE_NAMES)
def test_stencil_case(case_name: str) -> None:
    """Stencil sum over f=sin(lon)*cos(lat) matches Julia reference on a real grid."""
    rule = _load_rule()
    spec_case = _STENCIL_CASES[case_name]
    golden_case = _GOLDEN_STENCIL[case_name]

    grid_opts = spec_case["grid"]
    grid = grids.lat_lon(
        variant="regular",
        nlon=int(grid_opts["nlon"]),
        nlat=int(grid_opts["nlat"]),
        R=float(grid_opts["R"]),
        ghosts=0,
        dtype="float64",
        pole_policy="none",
        lon_start=0.0,
    )

    golden_by_ji = {
        (int(r["j"]), int(r["i"])): float(r["result"])
        for r in golden_case["results"]
    }

    for qp in spec_case["interior_query_points"]:
        j, i = int(qp[0]), int(qp[1])
        bindings = _bindings_for(grid, j, i)
        neighbors = grid.neighbors(j, i)

        stencil_sum = 0.0
        for entry in rule["stencil"]:
            sel = entry["selector"]
            coeff = eval_coeff(entry["coeff"], bindings)
            direction = _DIR_FOR_SELECTOR[(sel["axis"], int(sel["offset"]))]
            nbr = neighbors[direction]
            assert nbr is not None, (
                f"{case_name} ({j},{i}): {direction} neighbor is None (unexpected pole)"
            )
            nbr_lon, nbr_lat = grid.cell_center(nbr[0], nbr[1])
            stencil_sum += coeff * math.sin(nbr_lon) * math.cos(nbr_lat)

        expected = golden_by_ji[(j, i)]
        assert _close_rel(stencil_sum, expected), (
            f"{case_name} ({j},{i}): stencil sum {stencil_sum:.15g}, expected {expected:.15g}"
        )
