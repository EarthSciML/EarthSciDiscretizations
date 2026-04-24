"""Cross-language conformance for the arakawa grid family.

Verifies the Python accessor runtime produces accessor outputs that match
the committed Julia-reference golden within the tolerance declared in
``tests/conformance/grids/arakawa/fixtures.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_toolkit import grids
from earthsci_toolkit.grids.arakawa import CartesianBase

HARNESS_DIR = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "conformance"
    / "grids"
    / "arakawa"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
GOLDEN_DIR = HARNESS_DIR / "golden"
_SPEC = json.loads(FIXTURES_PATH.read_text())
REL_TOL = float(_SPEC["tolerance"]["relative"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]

_LOCATIONS = ("cell_center", "u_edge", "v_edge", "corner")
_STAGGERS = ("A", "B", "C", "D", "E")
_DIRS = ("W", "E", "S", "N")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _coord_at(g_a, g_c, loc: str, i: int, j: int) -> tuple[float, float]:
    # g_a is stagger=A (u/v at cell_center); g_c is stagger=C (u at u_edge,
    # v at v_edge). Use the stagger that routes u_face/v_face to the location
    # being queried, and cell_centers/corners for the other two.
    if loc == "cell_center":
        return g_a.cell_centers(i, j)
    if loc == "corner":
        return g_a.corners(i, j)
    if loc == "u_edge":
        return g_c.u_face(i, j)
    if loc == "v_edge":
        return g_c.v_face(i, j)
    raise AssertionError(f"unknown loc {loc!r}")


def _neighbors_wire(g, loc: str, i: int, j: int) -> dict:
    got = g.neighbors(loc, i, j)
    return {
        dkey: (None if got[dkey] is None else [int(got[dkey][0]), int(got[dkey][1])])
        for dkey in _DIRS
    }


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_arakawa_conformance_against_golden(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    opts = fixture["opts"]
    base_opts = opts["base"]
    assert base_opts["family"] == "cartesian", "only cartesian base supported"
    base = CartesianBase(
        xlo=float(base_opts["xlo"]),
        xhi=float(base_opts["xhi"]),
        ylo=float(base_opts["ylo"]),
        yhi=float(base_opts["yhi"]),
        nx=int(base_opts["nx"]),
        ny=int(base_opts["ny"]),
    )
    dtype = opts["dtype"]
    ghosts = int(opts["ghosts"])
    golden = json.loads((GOLDEN_DIR / f"{fixture_name}.json").read_text())

    assert base.nx * base.ny == golden["n_cells"]

    g_a = grids.arakawa(base=base, stagger="A", dtype=dtype, ghosts=ghosts)
    g_c = grids.arakawa(base=base, stagger="C", dtype=dtype, ghosts=ghosts)

    # --- coord + neighbor tables (stagger-independent) -------------------
    for lname in _LOCATIONS:
        ctable = golden["coords"][lname]
        ntable = golden["neighbors"][lname]
        for k, (i, j) in enumerate(ctable["points"]):
            x, y = _coord_at(g_a, g_c, lname, int(i), int(j))
            ex_x, ex_y = ctable["xy"][k]
            assert _close_rel(x, float(ex_x)), (
                f"coord {lname} x mismatch at qp[{k}]=({i},{j}): {x} vs {ex_x}"
            )
            assert _close_rel(y, float(ex_y)), (
                f"coord {lname} y mismatch at qp[{k}]=({i},{j}): {y} vs {ex_y}"
            )
            got = _neighbors_wire(g_a, lname, int(i), int(j))
            for dkey in _DIRS:
                expected = ntable[dkey][k]
                if expected is None:
                    assert got[dkey] is None, (
                        f"neighbor {dkey} {lname} at qp[{k}]=({i},{j}): "
                        f"expected None, got {got[dkey]}"
                    )
                else:
                    assert got[dkey] == [int(expected[0]), int(expected[1])], (
                        f"neighbor {dkey} {lname} at qp[{k}]=({i},{j}): "
                        f"{got[dkey]} vs {expected}"
                    )

    # --- metric_eval (stagger-independent for cartesian) -----------------
    mtable = golden["metrics"]
    for k, (i, j) in enumerate(mtable["points"]):
        for name in ("dx", "dy", "area"):
            got = g_a.metric_eval(name, int(i), int(j))
            expected = float(mtable[name][k])
            assert _close_rel(got, expected), (
                f"metric {name} mismatch at qp[{k}]=({i},{j}): {got} vs {expected}"
            )

    # --- per-stagger variable locations + shapes -------------------------
    for sname in fixture["staggers"]:
        assert sname in _STAGGERS
        g_s = grids.arakawa(base=base, stagger=sname, dtype=dtype, ghosts=ghosts)
        stab = golden["staggers"][sname]

        assert g_s.rotated == bool(stab["rotated"]), f"stagger={sname} rotated flag"
        for var in ("h", "u", "v"):
            got_loc = g_s.variable_location(var)
            assert got_loc == stab["variable_locations"][var], (
                f"stagger={sname} var={var} location"
            )
            got_shape = list(g_s.variable_shape(var))
            expected_shape = [int(x) for x in stab["variable_shapes"][var]]
            assert got_shape == expected_shape, (
                f"stagger={sname} var={var} shape"
            )
        for lname in _LOCATIONS:
            got_shape = list(g_s.location_shape(lname))
            expected_shape = [int(x) for x in stab["location_shapes"][lname]]
            assert got_shape == expected_shape, (
                f"stagger={sname} loc={lname} shape"
            )
