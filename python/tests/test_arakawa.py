"""Tests for the arakawa grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §7 contract:
stagger → variable-location table, location shapes, accessor coordinates
(C-grid recovers Cartesian offsets exactly), neighbor topology, metric
eval, and declarative .esm lowering shape.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import earthsci_toolkit
from earthsci_toolkit import grids
from earthsci_toolkit.grids import ArakawaGrid, BaseGrid, CartesianBase, arakawa
from earthsci_toolkit.grids.arakawa import (
    _LOCATIONS,
    _STAGGERS,
    _VARIABLE_LOCATIONS,
    _location_shape,
)

# ---------------------------------------------------------------- exports --


def test_exported_from_grids_namespace():
    assert grids.arakawa is arakawa
    assert grids.ArakawaGrid is ArakawaGrid
    assert grids.CartesianBase is CartesianBase


def test_base_grid_protocol_exported():
    assert isinstance(CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=2, ny=2), BaseGrid)


# ------------------------------------------- stagger → variable locations --


def test_stagger_variable_locations_match_spec():
    assert _VARIABLE_LOCATIONS["A"] == ("cell_center", "cell_center", "cell_center")
    assert _VARIABLE_LOCATIONS["B"] == ("cell_center", "corner", "corner")
    assert _VARIABLE_LOCATIONS["C"] == ("cell_center", "u_edge", "v_edge")
    assert _VARIABLE_LOCATIONS["D"] == ("cell_center", "v_edge", "u_edge")
    # E is topologically like B (45° rotation carried separately in .esm).
    assert _VARIABLE_LOCATIONS["E"] == ("cell_center", "corner", "corner")


def test_location_shape_table():
    assert _location_shape("cell_center", 10, 20) == (10, 20)
    assert _location_shape("u_edge", 10, 20) == (11, 20)
    assert _location_shape("v_edge", 10, 20) == (10, 21)
    assert _location_shape("corner", 10, 20) == (11, 21)


def test_all_staggers_covered():
    assert _STAGGERS == {"A", "B", "C", "D", "E"}
    for s in _STAGGERS:
        assert s in _VARIABLE_LOCATIONS


def test_all_locations_covered():
    assert _LOCATIONS == {"cell_center", "u_edge", "v_edge", "corner"}


# ------------------------------------------ CartesianBase validation --


def test_cartesian_base_rejects_zero_nx():
    with pytest.raises(ValueError):
        CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=0, ny=4)


def test_cartesian_base_rejects_negative_ny():
    with pytest.raises(ValueError):
        CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=4, ny=-1)


def test_cartesian_base_rejects_reversed_extent():
    with pytest.raises(ValueError):
        CartesianBase(xlo=1.0, xhi=0.0, ylo=0.0, yhi=1.0, nx=4, ny=4)
    with pytest.raises(ValueError):
        CartesianBase(xlo=0.0, xhi=1.0, ylo=1.0, yhi=1.0, nx=4, ny=4)


def test_cartesian_base_rejects_nonfinite():
    with pytest.raises(ValueError):
        CartesianBase(xlo=float("nan"), xhi=1.0, ylo=0.0, yhi=1.0, nx=4, ny=4)


def test_cartesian_base_rejects_non_int_nx():
    with pytest.raises(TypeError):
        CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=3.5, ny=4)  # type: ignore[arg-type]


def test_cartesian_base_basic_properties():
    b = CartesianBase(xlo=0.0, xhi=2.0, ylo=0.0, yhi=6.0, nx=4, ny=3)
    assert b.nx == 4 and b.ny == 3
    assert math.isclose(b.dx(), 0.5, rel_tol=0, abs_tol=1e-15)
    assert math.isclose(b.dy(), 2.0, rel_tol=0, abs_tol=1e-15)


# ------------------------------------------- generator / error contract --


def _cart(nx: int = 4, ny: int = 4) -> CartesianBase:
    return CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=nx, ny=ny)


def test_arakawa_keyword_only_signature():
    with pytest.raises(TypeError):
        arakawa(_cart(), "C")  # type: ignore[misc]


def test_arakawa_missing_base_raises_typeerror():
    with pytest.raises(TypeError):
        arakawa(stagger="C")  # type: ignore[call-arg]


def test_arakawa_missing_stagger_raises_typeerror():
    with pytest.raises(TypeError):
        arakawa(base=_cart())  # type: ignore[call-arg]


def test_arakawa_rejects_unknown_stagger():
    with pytest.raises(ValueError):
        arakawa(base=_cart(), stagger="Q")


def test_arakawa_rejects_non_base_object():
    with pytest.raises(TypeError):
        arakawa(base="not a grid", stagger="C")  # type: ignore[arg-type]


def test_arakawa_rejects_negative_ghosts():
    with pytest.raises(ValueError):
        arakawa(base=_cart(), stagger="C", ghosts=-1)


def test_arakawa_rejects_non_int_ghosts():
    with pytest.raises(TypeError):
        arakawa(base=_cart(), stagger="C", ghosts=1.5)  # type: ignore[arg-type]


def test_arakawa_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        arakawa(base=_cart(), stagger="C", dtype="float16")


def test_arakawa_defaults_match_contract():
    g = arakawa(base=_cart(), stagger="C")
    assert g.family == "arakawa"
    assert g.topology == "block_structured"
    assert g.dtype == "float64"
    assert g.ghosts == 0
    assert g.stagger == "C"
    assert g.nx == 4 and g.ny == 4
    assert g.n_cells == 16
    assert g.rotated is False


def test_arakawa_accepts_numpy_dtype():
    g = arakawa(base=_cart(), stagger="C", dtype=np.dtype("float32"))
    assert g.dtype == "float32"


# ---------------------------------- C-grid accessors recover cartesian --


def test_c_grid_accessors_recover_cartesian_coordinates():
    base = CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=10, ny=10)
    g = arakawa(base=base, stagger="C")

    cx, cy = g.cell_centers(0, 0)
    assert math.isclose(cx, 0.05, abs_tol=1e-14)
    assert math.isclose(cy, 0.05, abs_tol=1e-14)

    cx, cy = g.cell_centers(9, 9)
    assert math.isclose(cx, 0.95, abs_tol=1e-14)
    assert math.isclose(cy, 0.95, abs_tol=1e-14)

    # C: u lives at u-faces. u(0, 0) sits on western boundary at mid-y of row 0.
    ux, uy = g.u_face(0, 0)
    assert math.isclose(ux, 0.0, abs_tol=1e-14)
    assert math.isclose(uy, 0.05, abs_tol=1e-14)
    # u(10, 0) is the eastern boundary of the last column.
    ux, _ = g.u_face(10, 0)
    assert math.isclose(ux, 1.0, abs_tol=1e-14)

    # C: v lives at v-faces.
    vx, vy = g.v_face(0, 0)
    assert math.isclose(vx, 0.05, abs_tol=1e-14)
    assert math.isclose(vy, 0.0, abs_tol=1e-14)
    _, vy = g.v_face(0, 10)
    assert math.isclose(vy, 1.0, abs_tol=1e-14)

    assert g.corners(0, 0) == (0.0, 0.0)
    assert g.corners(10, 10) == (1.0, 1.0)


def test_c_grid_out_of_bounds_index_rejected():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.cell_centers(4, 0)  # nx=4
    with pytest.raises(ValueError):
        g.u_face(6, 0)  # u shape (5, 4) for nx=ny=4
    with pytest.raises(ValueError):
        g.v_face(0, 6)
    with pytest.raises(ValueError):
        g.corners(5, 5)  # corner shape (5, 5); valid is [0, 5)


# --------------------------------------- A-grid: colocation at centers --


def test_a_grid_colocates_u_v_at_centers():
    g = arakawa(base=_cart(), stagger="A")
    assert g.variable_shape("h") == (4, 4)
    assert g.variable_shape("u") == (4, 4)
    assert g.variable_shape("v") == (4, 4)
    assert g.u_face(1, 2) == g.cell_centers(1, 2)
    assert g.v_face(1, 2) == g.cell_centers(1, 2)


# --------------------------------------- B-grid: u, v at corners --


def test_b_grid_puts_uv_at_corners():
    g = arakawa(base=_cart(), stagger="B")
    assert g.variable_shape("h") == (4, 4)
    assert g.variable_shape("u") == (5, 5)
    assert g.variable_shape("v") == (5, 5)
    assert g.u_face(0, 0) == g.corners(0, 0)
    assert g.v_face(2, 1) == g.corners(2, 1)


# --------------------------------------- D-grid: swaps C's u/v faces --


def test_d_grid_swaps_c_uv_faces():
    base = CartesianBase(xlo=0.0, xhi=2.0, ylo=0.0, yhi=1.0, nx=4, ny=4)
    c = arakawa(base=base, stagger="C")
    d = arakawa(base=base, stagger="D")
    assert d.variable_shape("u") == c.variable_shape("v")
    assert d.variable_shape("v") == c.variable_shape("u")
    assert d.u_face(0, 0) == c.v_face(0, 0)
    assert d.v_face(0, 0) == c.u_face(0, 0)


# --------------------------------------- E-grid: topologically like B --


def test_e_grid_mirrors_b_topology_with_rotated_flag():
    b = arakawa(base=_cart(), stagger="B")
    e = arakawa(base=_cart(), stagger="E")
    assert e.variable_shape("h") == b.variable_shape("h")
    assert e.variable_shape("u") == b.variable_shape("u")
    assert e.variable_shape("v") == b.variable_shape("v")
    assert e.rotated is True
    assert b.rotated is False


# --------------------------------------- variable location lookup --


def test_variable_location_rejects_unknown_var():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.variable_location("w")


# --------------------------------------- neighbors --


def test_neighbors_interior_cell():
    g = arakawa(base=_cart(), stagger="C")
    n = g.cell_neighbors(1, 1)
    assert n["W"] == (0, 1)
    assert n["E"] == (2, 1)
    assert n["S"] == (1, 0)
    assert n["N"] == (1, 2)


def test_neighbors_corner_cell_has_none_on_boundary():
    g = arakawa(base=_cart(), stagger="C")
    n = g.cell_neighbors(0, 0)
    assert n["W"] is None
    assert n["S"] is None
    assert n["E"] == (1, 0)
    assert n["N"] == (0, 1)


def test_neighbors_bounds_check():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.cell_neighbors(4, 0)
    with pytest.raises(ValueError):
        g.cell_neighbors(0, 4)


def test_neighbors_uedge_last_column():
    g = arakawa(base=_cart(), stagger="C")
    n = g.neighbors("u_edge", 4, 3)  # shape (5, 4)
    assert n["E"] is None
    assert n["N"] is None
    assert n["W"] == (3, 3)
    assert n["S"] == (4, 2)


def test_neighbors_unknown_location_rejected():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.neighbors("face", 0, 0)


# --------------------------------------- coord accessor --


def test_coord_generic_accessor_matches_named_variants():
    g = arakawa(base=_cart(), stagger="C")
    assert g.coord("cell_center", 1, 2) == g.cell_centers(1, 2)
    assert g.coord("u_edge", 1, 2) == g.base.x_edge(1, 2)
    assert g.coord("v_edge", 1, 2) == g.base.y_edge(1, 2)
    assert g.coord("corner", 1, 2) == g.corners(1, 2)


def test_coord_rejects_unknown_location():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.coord("face", 0, 0)


# --------------------------------------- metric_eval --


def test_metric_eval_cartesian_uniform():
    base = CartesianBase(xlo=0.0, xhi=2.0, ylo=0.0, yhi=6.0, nx=4, ny=3)
    g = arakawa(base=base, stagger="C")
    assert math.isclose(g.metric_eval("dx", 0, 0), 0.5, abs_tol=1e-14)
    assert math.isclose(g.metric_eval("dy", 0, 0), 2.0, abs_tol=1e-14)
    assert math.isclose(g.metric_eval("area", 3, 2), 1.0, abs_tol=1e-14)


def test_metric_eval_rejects_unknown_name():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.metric_eval("not_a_metric", 0, 0)


def test_metric_eval_rejects_out_of_range_cell():
    g = arakawa(base=_cart(), stagger="C")
    with pytest.raises(ValueError):
        g.metric_eval("dx", 4, 0)


# --------------------------------------- to_esm: declarative config --


def test_to_esm_is_declarative():
    base = CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=5, ny=5)
    g = arakawa(base=base, stagger="C", ghosts=2)
    doc = g.to_esm()
    assert doc["family"] == "arakawa"
    assert doc["version"] == "1.0.0"
    assert doc["dtype"] == "float64"
    assert doc["topology"] == "block_structured"
    assert doc["stagger"] == "C"
    assert doc["ghosts"] == 2
    assert doc["n_cells"] == 25
    assert doc["rotated"] is False

    base_doc = doc["base"]
    assert base_doc["family"] == "cartesian"
    assert base_doc["nx"] == 5
    assert base_doc["ny"] == 5
    assert base_doc["extent"] == [[0.0, 0.0], [1.0, 1.0]]

    # Per §0 correction: no inline geometry arrays.
    import json

    wire = json.dumps(doc)
    assert '"cells":' not in wire
    assert '"edges":' not in wire
    assert '"vertices":' not in wire
    assert len(wire) < 2_000


def test_to_esm_rotated_flag_for_e():
    base = CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=4, ny=4)
    e = arakawa(base=base, stagger="E")
    b = arakawa(base=base, stagger="B")
    assert e.to_esm()["rotated"] is True
    assert b.to_esm()["rotated"] is False


def test_provenance_identifies_python_binding():
    g = arakawa(base=_cart(), stagger="C")
    prov = g.to_esm()["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_toolkit.__version__
    assert prov["stagger"] == "C"
    assert prov["source"] == "earthsci_toolkit.grids.arakawa"


def test_dtype_f32_propagates_to_esm():
    g = arakawa(base=_cart(), stagger="C", dtype="float32")
    assert g.dtype == "float32"
    assert g.to_esm()["dtype"] == "float32"


# --------------------------------------- cross-stagger shape consistency --


@pytest.mark.parametrize("stagger", ["A", "B", "C", "D", "E"])
def test_variable_shapes_consistent_with_stagger_table(stagger):
    base = CartesianBase(xlo=0.0, xhi=1.0, ylo=0.0, yhi=1.0, nx=7, ny=5)
    g = arakawa(base=base, stagger=stagger)
    h_loc, u_loc, v_loc = _VARIABLE_LOCATIONS[stagger]
    assert g.variable_shape("h") == _location_shape(h_loc, 7, 5)
    assert g.variable_shape("u") == _location_shape(u_loc, 7, 5)
    assert g.variable_shape("v") == _location_shape(v_loc, 7, 5)


# --------------------------------------- repr --


def test_repr_contains_stagger_and_shape():
    g = arakawa(base=_cart(), stagger="C")
    r = repr(g)
    assert "ArakawaGrid" in r
    assert "'C'" in r
    assert "nx=4" in r
    assert "ny=4" in r
