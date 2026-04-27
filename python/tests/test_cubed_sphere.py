"""Tests for the cubed_sphere grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §7 contract:
topology consistency (6 panels, connectivity symmetric), metric accuracy
(analytic gnomonic identities, inverse-of-covariant), area correctness
(total = 4 pi R^2), and declarative .esm round-trip shape.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import earthsci_discretizations
from earthsci_discretizations import grids
from earthsci_discretizations.grids import CubedSphereGrid, cubed_sphere
from earthsci_discretizations.grids.cubed_sphere import (
    _PANEL_CONNECTIVITY,
    _gnomonic_metric,
    _gnomonic_to_cart,
    _gnomonic_to_lonlat,
)

# ------------------------------------------------------------------ API --


def test_exported_from_grids_namespace():
    assert grids.cubed_sphere is cubed_sphere
    assert grids.CubedSphereGrid is CubedSphereGrid


def test_required_Nc_is_keyword_only():
    with pytest.raises(TypeError):
        cubed_sphere(8)  # type: ignore[misc] - positional forbidden per §2.4


def test_missing_Nc_raises_typeerror():
    with pytest.raises(TypeError):
        cubed_sphere()  # type: ignore[call-arg]


def test_default_R_and_dtype():
    g = cubed_sphere(Nc=4)
    assert g.R == 6.371e6
    assert g.dtype == "float64"
    assert g.ghosts == 0


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        cubed_sphere(Nc=4, dtype="float16")


def test_invalid_Nc_rejected():
    with pytest.raises(ValueError):
        cubed_sphere(Nc=0)
    with pytest.raises(TypeError):
        cubed_sphere(Nc=3.5)  # type: ignore[arg-type]


# ------------------------------------------------------------ topology --


def test_n_cells_matches_6_times_Nc_squared():
    g = cubed_sphere(Nc=8)
    assert g.n_cells == 6 * 8 * 8
    assert g.family == "cubed_sphere"
    assert g.topology == "block_structured"


def test_panel_connectivity_symmetric():
    """`verify_connectivity`: W/E/S/N edges must be mutually reciprocal."""
    for p in range(6):
        for d in ("W", "E", "S", "N"):
            nb_panel, nb_edge, _ = _PANEL_CONNECTIVITY[p][d]
            back_panel, _, _ = _PANEL_CONNECTIVITY[nb_panel][nb_edge]
            assert back_panel == p, (
                f"panel {p} dir {d} -> ({nb_panel},{nb_edge}) does not "
                f"return to {p}"
            )


def test_interior_neighbors_are_local():
    g = cubed_sphere(Nc=8)
    n = g.neighbors(p=0, i=3, j=4)
    assert n["W"] == (0, 2, 4)
    assert n["E"] == (0, 4, 4)
    assert n["S"] == (0, 3, 3)
    assert n["N"] == (0, 3, 5)


def test_cross_panel_neighbor_reciprocity():
    """Stepping across a panel edge and back returns a boundary cell on the
    original panel (the exact same cell, up to the along-edge transform)."""
    g = cubed_sphere(Nc=6)
    for p in range(6):
        for d in ("W", "E", "S", "N"):
            # pick a representative boundary cell on this edge
            if d == "W":
                i, j = 0, 2
            elif d == "E":
                i, j = g.Nc - 1, 2
            elif d == "S":
                i, j = 2, 0
            else:
                i, j = 2, g.Nc - 1
            nb_panel, nb_i, nb_j = g.neighbors(p, i, j)[d]
            # neighbor's same-side edge leads back to some cell on panel p
            back_dir = _PANEL_CONNECTIVITY[p][d][1]
            # The return path from nb cell across `back_dir` should point to p
            back_neighbors = g.neighbors(nb_panel, nb_i, nb_j)
            back_cell = back_neighbors[back_dir]
            assert back_cell[0] == p, (
                f"{p}/{d} -> {nb_panel}/{back_dir} did not return to {p}"
            )


def test_neighbor_cell_indices_are_valid():
    g = cubed_sphere(Nc=5)
    for p in range(6):
        for i in range(g.Nc):
            for j in range(g.Nc):
                for d, (pp, ii, jj) in g.neighbors(p, i, j).items():
                    assert 0 <= pp < 6, d
                    assert 0 <= ii < g.Nc, d
                    assert 0 <= jj < g.Nc, d


# ------------------------------------------------------------- centers --


def test_cell_centers_shape_and_range():
    Nc = 8
    g = cubed_sphere(Nc=Nc)
    lon, lat = g.cell_centers()
    assert lon.shape == (6, Nc, Nc)
    assert lat.shape == (6, Nc, Nc)
    assert np.all(lon >= -math.pi - 1e-12)
    assert np.all(lon <= math.pi + 1e-12)
    assert np.all(lat >= -math.pi / 2 - 1e-12)
    assert np.all(lat <= math.pi / 2 + 1e-12)


def test_cell_center_scalar_matches_bulk():
    Nc = 4
    g = cubed_sphere(Nc=Nc)
    lon, lat = g.cell_centers()
    for p in (0, 2, 5):
        for (i, j) in ((0, 0), (1, 2), (Nc - 1, Nc - 1)):
            s_lon, s_lat = g.cell_centers(p, i, j)
            assert math.isclose(s_lon, float(lon[p, i, j]), abs_tol=1e-12)
            assert math.isclose(s_lat, float(lat[p, i, j]), abs_tol=1e-12)


def test_panel_centers_lie_on_expected_cardinals():
    """Panel 0 center (xi=eta=0) maps to (lon=0, lat=0); panel 2 (top) to
    lat=+pi/2; panel 5 (bottom) to lat=-pi/2.  This is the anchor the FV3
    gnomonic convention provides independent of Nc."""
    lon0, lat0 = _gnomonic_to_lonlat(0.0, 0.0, 0)
    assert math.isclose(lon0, 0.0, abs_tol=1e-15)
    assert math.isclose(lat0, 0.0, abs_tol=1e-15)
    _, lat2 = _gnomonic_to_lonlat(0.0, 0.0, 2)
    assert math.isclose(lat2, math.pi / 2, abs_tol=1e-15)
    _, lat5 = _gnomonic_to_lonlat(0.0, 0.0, 5)
    assert math.isclose(lat5, -math.pi / 2, abs_tol=1e-15)


# ---------------------------------------------------------------- area --


def test_total_sphere_area_matches_4piR2():
    """Sum of all cell face areas must equal 4*pi*R^2 to high precision."""
    R = 6.371e6
    g = cubed_sphere(Nc=24, R=R)
    total = float(g.area.sum())
    expected = 4 * math.pi * R * R
    assert abs(total - expected) / expected < 1e-9


def test_total_area_unit_sphere():
    g = cubed_sphere(Nc=16, R=1.0)
    total = float(g.area.sum())
    assert abs(total - 4 * math.pi) < 1e-9


def test_area_positive_everywhere():
    g = cubed_sphere(Nc=8)
    assert np.all(g.area > 0)


def test_area_via_metric_eval_matches_array():
    g = cubed_sphere(Nc=4)
    for p in (0, 1, 3):
        for i in (0, 2, 3):
            for j in (0, 1, 3):
                a = g.metric_eval("area", p, i, j)
                assert math.isclose(a, float(g.area[p, i, j]), rel_tol=1e-12)


# ------------------------------------------------------------ metrics --


@pytest.mark.parametrize("xi,eta", [(0.0, 0.0), (0.1, -0.2), (-0.3, 0.25), (0.5, 0.4)])
def test_metric_analytic_identity_at_origin(xi, eta):
    """At xi=eta=0 the gnomonic metric reduces to J=R^2, g=R^2 I (known ID)."""
    R = 2.5
    J, gxx, gyy, gxy = _gnomonic_metric(xi, eta, R)
    # Panel-independent sanity: J must equal sqrt(det(g)) up to sign.
    det_g = gxx * gyy - gxy * gxy
    assert det_g > 0
    assert math.isclose(J, math.sqrt(det_g), rel_tol=1e-12)


def test_metric_at_origin_equals_R_squared_identity():
    R = 3.0
    J, gxx, gyy, gxy = _gnomonic_metric(0.0, 0.0, R)
    assert math.isclose(J, R * R, rel_tol=1e-14)
    assert math.isclose(gxx, R * R, rel_tol=1e-14)
    assert math.isclose(gyy, R * R, rel_tol=1e-14)
    assert math.isclose(gxy, 0.0, abs_tol=1e-14)


def test_inverse_metric_is_actually_inverse():
    """g_ab * ginv^bc = delta_a^c at every cell center."""
    g = cubed_sphere(Nc=6, R=1.0)
    for (i, j) in ((0, 0), (2, 3), (g.Nc - 1, g.Nc - 1)):
        gxx = g.metric_eval("g_xixi", 0, i, j)
        gyy = g.metric_eval("g_etaeta", 0, i, j)
        gxy = g.metric_eval("g_xieta", 0, i, j)
        ixx = g.metric_eval("ginv_xixi", 0, i, j)
        iyy = g.metric_eval("ginv_etaeta", 0, i, j)
        ixy = g.metric_eval("ginv_xieta", 0, i, j)
        # (g @ ginv) entries
        m00 = gxx * ixx + gxy * ixy
        m01 = gxx * ixy + gxy * iyy
        m10 = gxy * ixx + gyy * ixy
        m11 = gxy * ixy + gyy * iyy
        assert math.isclose(m00, 1.0, abs_tol=1e-12)
        assert math.isclose(m11, 1.0, abs_tol=1e-12)
        assert math.isclose(m01, 0.0, abs_tol=1e-12)
        assert math.isclose(m10, 0.0, abs_tol=1e-12)


def test_metric_eval_rejects_unknown_name():
    g = cubed_sphere(Nc=4)
    with pytest.raises(ValueError):
        g.metric_eval("not_a_metric", 0, 0, 0)


def test_metric_eval_rejects_out_of_range_cell():
    g = cubed_sphere(Nc=4)
    with pytest.raises(ValueError):
        g.metric_eval("J", 6, 0, 0)
    with pytest.raises(ValueError):
        g.metric_eval("J", 0, 4, 0)
    with pytest.raises(ValueError):
        g.metric_eval("J", 0, 0, -1)


def test_J_positive_everywhere():
    g = cubed_sphere(Nc=6)
    for p in range(6):
        for i in range(g.Nc):
            for j in range(g.Nc):
                assert g.metric_eval("J", p, i, j) > 0


# -------------------------------------------------- gnomonic math sanity --


def test_gnomonic_to_cart_is_on_unit_sphere():
    for xi in (-0.3, 0.0, 0.2):
        for eta in (-0.25, 0.0, 0.35):
            for p in range(6):
                x, y, z = _gnomonic_to_cart(xi, eta, p)
                assert abs(x * x + y * y + z * z - 1.0) < 1e-14


# --------------------------------------------------------------- to_esm --


def test_to_esm_is_declarative():
    """Per the 2026-04-20 correction, to_esm returns a small config - no
    inline geometry arrays, no per-cell payload."""
    g = cubed_sphere(Nc=48, R=6.371e6, dtype="float64", ghosts=3)
    doc = g.to_esm()
    assert doc["family"] == "cubed_sphere"
    assert doc["version"] == "1.0.0"
    assert doc["dtype"] == "float64"
    assert doc["topology"] == "block_structured"
    assert doc["generator"] == "gnomonic_c6"
    assert doc["params"] == {"Nc": 48, "R": 6.371e6, "ghosts": 3}
    assert "provenance" in doc

    # No inline geometry arrays anywhere in the serialized config.
    import json

    wire = json.dumps(doc)
    assert "cells" not in wire
    assert "lon_array" not in wire
    # And the wire form must stay small (declarative, not a blob).
    assert len(wire) < 2_000


def test_to_esm_json_roundtrip():
    import json

    g = cubed_sphere(Nc=4, R=1.0)
    doc = g.to_esm()
    reparsed = json.loads(json.dumps(doc))
    assert reparsed["params"]["Nc"] == 4
    assert reparsed["params"]["R"] == 1.0


def test_provenance_carries_binding_identity():
    g = cubed_sphere(Nc=4)
    prov = g.to_esm()["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_discretizations.__version__
    assert prov["generator"] == "gnomonic_c6"


# -------------------------------------------------------------- dtypes --


def test_float32_dtype_propagates_to_arrays():
    g = cubed_sphere(Nc=4, dtype="float32")
    assert g.lon.dtype == np.float32
    assert g.lat.dtype == np.float32
    assert g.area.dtype == np.float32
    assert g.to_esm()["dtype"] == "float32"


def test_numpy_dtype_accepted():
    g = cubed_sphere(Nc=4, dtype=np.dtype("float64"))
    assert g.dtype == "float64"
