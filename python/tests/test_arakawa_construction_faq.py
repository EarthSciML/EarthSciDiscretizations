"""Cross-language conformance for the arakawa **construction FAQ** (esd-dru).

Builds each fixture's cartesian base, re-derives the staggered construction
through the elementwise FAQ bridge (:func:`arakawa_construction_faq`, routed via
the ESS ``eval_coeff`` evaluator), and asserts the result is

1. byte-identical to the imperative ``CartesianBase`` coordinate accessors over
   **all** cells of every location (parity — the FAQ introduces zero error), and
2. equal to the committed Julia-reference
   ``tests/conformance/grids/arakawa/construction/golden.json`` (the cartesian
   base geometry is pure affine arithmetic, well within the declared 1e-14).

The construction golden samples the four staggered location coordinate arrays,
the cell-center neighbor maps, and the boundary masks at curated ``points``, and
carries the full static A–E variable-location / shape table (stagger-independent
geometry; only the variable-location mapping depends on the stagger).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_discretizations.grids.arakawa import CartesianBase
from earthsci_discretizations.grids.arakawa_faq import arakawa_construction_faq

HARNESS_DIR = (
    Path(__file__).resolve().parents[2] / "tests" / "conformance" / "grids" / "arakawa"
)
_SPEC = json.loads((HARNESS_DIR / "fixtures.json").read_text())
_GOLDEN = json.loads((HARNESS_DIR / "construction" / "golden.json").read_text())
REL_TOL = float(_GOLDEN["tolerance"]["relative"])
_BY_NAME = {f["name"]: f for f in _GOLDEN["fixtures"]}
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]
_LOCATIONS = ("cell_center", "u_edge", "v_edge", "corner")
_STAGGERS = ("A", "B", "C", "D", "E")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _loc_shape(loc: str, nx: int, ny: int) -> tuple[int, int]:
    return {
        "cell_center": (nx, ny),
        "u_edge": (nx + 1, ny),
        "v_edge": (nx, ny + 1),
        "corner": (nx + 1, ny + 1),
    }[loc]


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_arakawa_construction_faq(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    bo = fixture["opts"]["base"]
    assert bo["family"] == "cartesian"
    base = CartesianBase(
        xlo=float(bo["xlo"]),
        xhi=float(bo["xhi"]),
        ylo=float(bo["ylo"]),
        yhi=float(bo["yhi"]),
        nx=int(bo["nx"]),
        ny=int(bo["ny"]),
    )
    gl = _BY_NAME[fixture_name]
    nx, ny = base.nx, base.ny
    faq = arakawa_construction_faq(base, "A")  # geometry is stagger-independent

    # --- counts + scalars ---
    assert faq.nx == int(gl["nx"]) and faq.ny == int(gl["ny"])
    assert nx * ny == int(gl["n_cells"])
    assert _close_rel(faq.dx, float(gl["dx"]))
    assert _close_rel(faq.dy, float(gl["dy"]))
    assert _close_rel(faq.cell_volume, float(gl["cell_volume"]))

    coords_faq = {
        "cell_center": faq.coords_cell_center,
        "u_edge": faq.coords_u_edge,
        "v_edge": faq.coords_v_edge,
        "corner": faq.coords_corner,
    }
    base_acc = {
        "cell_center": base.cell_center,
        "u_edge": base.x_edge,
        "v_edge": base.y_edge,
        "corner": base.corner,
    }

    # --- parity: FAQ == imperative base accessors, over ALL cells/locations ---
    for loc in _LOCATIONS:
        ni, nj = _loc_shape(loc, nx, ny)
        xs, ys = coords_faq[loc]
        acc = base_acc[loc]
        for j in range(nj):
            for i in range(ni):
                ex, ey = acc(i, j)
                k = j * ni + i
                assert float(xs[k]) == float(ex), f"{fixture_name}:{loc} x parity ({i},{j})"
                assert float(ys[k]) == float(ey), f"{fixture_name}:{loc} y parity ({i},{j})"

    # --- cross-binding golden: coords at curated points ---
    for loc in _LOCATIONS:
        ni, _ = _loc_shape(loc, nx, ny)
        pts = gl["coords"]["points"][loc]
        gxy = json.loads(gl["coords"][loc])
        xs, ys = coords_faq[loc]
        for k, (i, j) in enumerate(pts):
            flat = int(j) * ni + int(i)
            assert _close_rel(float(xs[flat]), float(gxy[k][0])), f"{fixture_name}:{loc} x golden"
            assert _close_rel(float(ys[flat]), float(gxy[k][1])), f"{fixture_name}:{loc} y golden"

    # --- cross-binding golden: cell-center boundary masks at curated points ---
    bd = gl["boundary"]
    bd_faq = {
        "x_lower": faq.boundary_x_lower,
        "x_upper": faq.boundary_x_upper,
        "y_lower": faq.boundary_y_lower,
        "y_upper": faq.boundary_y_upper,
    }
    for key, arr in bd_faq.items():
        exp = json.loads(bd[key])
        for k, (i, j) in enumerate(bd["points"]):
            assert int(arr[int(j) * nx + int(i)]) == int(exp[k]), f"{fixture_name}: boundary {key}"

    # --- cross-binding golden: cell-center neighbor maps at curated points ---
    ni_g = gl["neighbor_indices"]
    nb_faq = {
        "x_minus": faq.neighbor_x_minus,
        "x_plus": faq.neighbor_x_plus,
        "y_minus": faq.neighbor_y_minus,
        "y_plus": faq.neighbor_y_plus,
    }
    for key, arr in nb_faq.items():
        exp = json.loads(ni_g[key])
        for k, (i, j) in enumerate(ni_g["points"]):
            assert int(arr[int(j) * nx + int(i)]) == int(exp[k]), f"{fixture_name}: neighbor {key}"

    # --- cross-binding golden: static A–E variable-location / shape table ---
    loc_shapes = {loc: list(_loc_shape(loc, nx, ny)) for loc in _LOCATIONS}
    for s in _STAGGERS:
        faq_s = arakawa_construction_faq(base, s)
        gs = gl["staggers"][s]
        assert faq_s.rotated == bool(gs["rotated"]), f"{fixture_name}: {s} rotated"
        h, u, v = faq_s.variable_locations
        assert {"h": h, "u": u, "v": v} == gs["variable_locations"], f"{fixture_name}: {s} var_loc"
        for loc in _LOCATIONS:
            assert loc_shapes[loc] == list(gs["location_shapes"][loc]), (
                f"{fixture_name}: {s} location_shapes {loc}"
            )
        for var, vloc in (("h", h), ("u", u), ("v", v)):
            assert loc_shapes[vloc] == list(gs["variable_shapes"][var]), (
                f"{fixture_name}: {s} variable_shapes {var}"
            )
