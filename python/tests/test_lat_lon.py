"""Tests for the lat_lon grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §7 contract:
topology consistency (periodic lon, polar boundaries under ``pole_policy='none'``),
metric accuracy (analytic identities, inverse-of-covariant), area correctness
(total = 4 pi R^2), declarative .esm round-trip shape, and reduced-Gaussian
variant support with user-supplied ``nlon_per_row`` / ``lat_edges``.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

import earthsci_discretizations
from earthsci_discretizations import grids
from earthsci_discretizations.grids import LatLonGrid, lat_lon
from earthsci_discretizations.grids.lat_lon import _map_i

# ------------------------------------------------------------------ API --


def test_exported_from_grids_namespace():
    assert grids.lat_lon is lat_lon
    assert grids.LatLonGrid is LatLonGrid


def test_required_nlon_is_keyword_only():
    with pytest.raises(TypeError):
        lat_lon(8, 4)  # type: ignore[misc] - positional forbidden per §2.4


def test_missing_nlon_raises_typeerror():
    with pytest.raises(TypeError):
        lat_lon(nlat=4)


def test_missing_nlat_raises_typeerror():
    with pytest.raises(TypeError):
        lat_lon(nlon=4)


def test_defaults_match_contract():
    g = lat_lon(nlon=8, nlat=4)
    assert g.variant == "regular"
    assert g.R == 6.371e6
    assert g.dtype == "float64"
    assert g.ghosts == 0
    assert g.family == "lat_lon"
    assert g.topology == "rectilinear"
    assert g.pole_policy == "none"
    assert g.nlon_uniform == 8
    assert g.nlat == 4
    assert g.n_cells == 32
    assert g.lon_start == pytest.approx(-math.pi)


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, dtype="float16")


def test_invalid_nlon_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=0, nlat=4)
    with pytest.raises(TypeError):
        lat_lon(nlon=3.5, nlat=4)  # type: ignore[arg-type]


def test_invalid_nlat_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=0)


def test_invalid_r_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, R=0)
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, R=-1.0)
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, R=float("nan"))


def test_invalid_variant_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, variant="bogus")


def test_invalid_pole_policy_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=4, pole_policy="spin")


def test_non_none_pole_policy_declared_but_not_implemented():
    with pytest.raises(ValueError, match="not implemented"):
        lat_lon(nlon=4, nlat=4, pole_policy="average")
    with pytest.raises(ValueError, match="not implemented"):
        lat_lon(nlon=4, nlat=4, pole_policy="fold")


def test_regular_rejects_nlon_per_row():
    with pytest.raises(ValueError):
        lat_lon(nlon=4, nlat=2, nlon_per_row=[4, 4])


def test_reduced_gaussian_rejects_nlon():
    with pytest.raises(ValueError):
        lat_lon(variant="reduced_gaussian", nlon=4, nlon_per_row=[4, 8])


def test_reduced_gaussian_requires_nlon_per_row():
    with pytest.raises(TypeError):
        lat_lon(variant="reduced_gaussian", nlat=4)


def test_reduced_gaussian_length_mismatch():
    with pytest.raises(ValueError):
        lat_lon(variant="reduced_gaussian", nlat=3, nlon_per_row=[4, 8])


def test_reduced_gaussian_rejects_zero_row():
    with pytest.raises(ValueError):
        lat_lon(variant="reduced_gaussian", nlon_per_row=[4, 0, 4])


def test_invalid_lat_edges_rejected():
    with pytest.raises(ValueError):
        lat_lon(nlon=2, nlat=2, lat_edges=[0.0, 0.0, 0.5])
    with pytest.raises(ValueError):
        lat_lon(nlon=2, nlat=1, lat_edges=[-2.0, 2.0])


def test_lat_centers_outside_edges_rejected():
    with pytest.raises(ValueError):
        lat_lon(
            nlon=2,
            nlat=2,
            lat_edges=[-0.5, 0.0, 0.5],
            lat_centers=[1.0, 0.25],
        )


# -------------------------------------------------------------- topology --


def test_lat_edges_default_span_pole_to_pole():
    g = lat_lon(nlon=4, nlat=6)
    assert g.lat_edges.shape == (7,)
    assert g.lat_edges[0] == pytest.approx(-math.pi / 2)
    assert g.lat_edges[-1] == pytest.approx(math.pi / 2)


def test_interior_neighbors_local():
    g = lat_lon(nlon=8, nlat=6)
    ns = g.neighbors(2, 3)
    assert ns["W"] == (2, 2)
    assert ns["E"] == (2, 4)
    assert ns["S"] == (1, 3)
    assert ns["N"] == (3, 3)


def test_longitude_wraps_periodically():
    g = lat_lon(nlon=8, nlat=4)
    assert g.neighbors(2, 0)["W"] == (2, 7)
    assert g.neighbors(2, 7)["E"] == (2, 0)


def test_south_pole_has_no_s_neighbor():
    g = lat_lon(nlon=8, nlat=4)
    ns = g.neighbors(0, 3)
    assert ns["S"] is None
    assert ns["N"] == (1, 3)


def test_north_pole_has_no_n_neighbor():
    g = lat_lon(nlon=8, nlat=4)
    ns = g.neighbors(3, 3)
    assert ns["N"] is None
    assert ns["S"] == (2, 3)


def test_all_neighbors_in_range():
    g = lat_lon(nlon=5, nlat=4)
    for j in range(g.nlat):
        for i in range(g.nlon(j)):
            for cell in g.neighbors(j, i).values():
                if cell is None:
                    continue
                jj, ii = cell
                assert 0 <= jj < g.nlat
                assert 0 <= ii < g.nlon(jj)


def test_reduced_gaussian_neighbor_rounding():
    g = lat_lon(variant="reduced_gaussian", nlon_per_row=[4, 8, 4])
    # Row 1 cell (1, 5) -> center frac = (5 + 0.5) / 8 = 0.6875
    # -> row 0 (width 4) column floor(0.6875 * 4) = 2.
    assert g.neighbors(1, 5)["S"] == (0, 2)
    assert g.neighbors(1, 0)["S"] == (0, 0)


def test_reduced_gaussian_cell_counts():
    g = lat_lon(variant="reduced_gaussian", nlon_per_row=[4, 8, 4])
    assert g.n_cells == 16
    assert g.nlat == 3
    assert g.row_offset(0) == 0
    assert g.row_offset(1) == 4
    assert g.row_offset(2) == 12


def test_reduced_gaussian_nlon_uniform_is_none():
    g = lat_lon(variant="reduced_gaussian", nlon_per_row=[4, 8, 4])
    assert g.nlon_uniform is None


def test_out_of_range_cell_errors():
    g = lat_lon(nlon=4, nlat=4)
    with pytest.raises(ValueError):
        g.neighbors(4, 0)
    with pytest.raises(ValueError):
        g.neighbors(0, 4)
    with pytest.raises(ValueError):
        g.cell_center(-1, 0)


# --------------------------------------------------------------- centers --


def test_cell_center_at_equator_equals_row_center_lat():
    # nlat=2 -> edges at [-pi/2, 0, pi/2]; row 0 center = -pi/4, row 1 = pi/4.
    g = lat_lon(nlon=4, nlat=2)
    _, lat0 = g.cell_center(0, 0)
    _, lat1 = g.cell_center(1, 0)
    assert lat0 == pytest.approx(-math.pi / 4)
    assert lat1 == pytest.approx(math.pi / 4)


def test_cell_centers_bulk_matches_scalar():
    g = lat_lon(nlon=6, nlat=4)
    lon_bulk, lat_bulk = g.cell_centers()
    for j in (0, 1, 3):
        for i in (0, 2, 5):
            s_lon, s_lat = g.cell_center(j, i)
            idx = g.row_offset(j) + i
            assert lon_bulk[idx] == pytest.approx(s_lon)
            assert lat_bulk[idx] == pytest.approx(s_lat)


def test_cell_centers_lon_lat_in_range():
    g = lat_lon(nlon=12, nlat=6)
    lon_bulk, lat_bulk = g.cell_centers()
    assert np.all(lon_bulk >= -math.pi - 1e-12)
    assert np.all(lon_bulk <= math.pi + 1e-12)
    assert np.all(lat_bulk >= -math.pi / 2 - 1e-12)
    assert np.all(lat_bulk <= math.pi / 2 + 1e-12)


def test_lon_lat_properties_are_cached():
    g = lat_lon(nlon=4, nlat=3)
    assert g.lon is g.lon
    assert g.lat is g.lat


def test_lon_edges_and_centers_lengths():
    g = lat_lon(nlon=5, nlat=3)
    assert g.lon_edges(0).shape == (6,)
    assert g.lon_centers(0).shape == (5,)


# ------------------------------------------------------------------ area --


def test_total_sphere_area_matches_4_pi_r_squared():
    r = 6.371e6
    g = lat_lon(nlon=72, nlat=36, R=r)
    total = float(g.area.sum())
    expected = 4.0 * math.pi * r * r
    assert abs(total - expected) / expected < 1e-12


def test_total_unit_sphere_area_is_4pi():
    g = lat_lon(nlon=10, nlat=8, R=1.0)
    total = float(g.area.sum())
    assert abs(total - 4.0 * math.pi) < 1e-12


def test_reduced_gaussian_total_area_is_4pi_r_squared():
    g = lat_lon(
        variant="reduced_gaussian",
        nlon_per_row=[4, 8, 12, 12, 8, 4],
        R=1.0,
    )
    total = float(g.area.sum())
    assert abs(total - 4.0 * math.pi) < 1e-12


def test_area_positive_everywhere():
    g = lat_lon(nlon=6, nlat=4)
    assert np.all(g.area > 0.0)


def test_area_via_metric_matches_bulk():
    g = lat_lon(nlon=6, nlat=4)
    bulk = g.area
    for j in (0, 1, 3):
        for i in (0, 2, 5):
            m = g.metric_eval("area", j, i)
            idx = g.row_offset(j) + i
            assert m == pytest.approx(bulk[idx], rel=1e-12)


# ---------------------------------------------------------------- metric --


def test_metric_at_equator_row_is_r_squared():
    # Many rows push the equator row latitude arbitrarily close to 0.
    r = 3.0
    g = lat_lon(nlon=4, nlat=400, R=r)
    j_eq = 200
    g_ll = g.metric_eval("g_lonlon", j_eq, 0)
    g_pp = g.metric_eval("g_latlat", j_eq, 0)
    assert abs(g_ll - r * r) / (r * r) < 1e-4
    assert g_pp == pytest.approx(r * r)


def test_metric_g_lonlat_is_zero():
    g = lat_lon(nlon=4, nlat=4)
    for j in range(g.nlat):
        assert g.metric_eval("g_lonlat", j, 0) == 0.0
        assert g.metric_eval("ginv_lonlat", j, 0) == 0.0


def test_inverse_metric_is_inverse():
    g = lat_lon(nlon=6, nlat=6, R=1.0)
    for j in (0, 2, 5):
        g_ll = g.metric_eval("g_lonlon", j, 0)
        g_pp = g.metric_eval("g_latlat", j, 0)
        i_ll = g.metric_eval("ginv_lonlon", j, 0)
        i_pp = g.metric_eval("ginv_latlat", j, 0)
        assert g_ll * i_ll == pytest.approx(1.0, rel=1e-12)
        assert g_pp * i_pp == pytest.approx(1.0, rel=1e-12)


def test_metric_eval_rejects_out_of_range_cell():
    g = lat_lon(nlon=4, nlat=4)
    with pytest.raises(ValueError):
        g.metric_eval("J", 4, 0)
    with pytest.raises(ValueError):
        g.metric_eval("J", 0, 4)


def test_metric_eval_rejects_unknown_name():
    g = lat_lon(nlon=4, nlat=4)
    with pytest.raises(ValueError, match="unknown metric"):
        g.metric_eval("not_a_metric", 0, 0)


def test_jacobian_positive_everywhere():
    g = lat_lon(nlon=6, nlat=6)
    for j in range(g.nlat):
        for i in range(g.nlon(j)):
            assert g.metric_eval("J", j, i) > 0.0


# ------------------------------------------------------------------ esm --


def test_to_esm_regular_is_declarative():
    g = lat_lon(nlon=72, nlat=36, R=6.371e6, ghosts=2)
    doc = g.to_esm()
    assert doc["family"] == "lat_lon"
    assert doc["version"] == "1.0.0"
    assert doc["dtype"] == "float64"
    assert doc["topology"] == "rectilinear"
    assert doc["variant"] == "regular"
    assert doc["generator"] == "lat_lon_regular"
    assert doc["params"]["nlon"] == 72
    assert doc["params"]["nlat"] == 36
    assert doc["params"]["R"] == 6.371e6
    assert doc["params"]["ghosts"] == 2
    assert doc["params"]["pole_policy"] == "none"
    assert "provenance" in doc

    wire = json.dumps(doc)
    assert "cells" not in wire
    assert "lon_array" not in wire
    assert len(wire) < 2_000


def test_to_esm_reduced_gaussian_carries_schedule():
    g = lat_lon(variant="reduced_gaussian", nlon_per_row=[4, 8, 4], R=1.0)
    doc = g.to_esm()
    assert doc["variant"] == "reduced_gaussian"
    assert doc["generator"] == "lat_lon_reduced_gaussian"
    assert doc["params"]["nlat"] == 3
    assert doc["params"]["nlon_per_row"] == [4, 8, 4]
    assert isinstance(doc["params"]["lat_edges"], list)


def test_to_esm_roundtrips_via_json():
    g = lat_lon(nlon=4, nlat=3, R=1.0)
    doc = g.to_esm()
    text = json.dumps(doc)
    reparsed = json.loads(text)
    assert reparsed["params"]["nlon"] == 4
    assert reparsed["params"]["nlat"] == 3
    assert reparsed["params"]["R"] == 1.0


def test_provenance_identifies_python_binding():
    g = lat_lon(nlon=4, nlat=4)
    prov = g.to_esm()["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_discretizations.__version__
    assert prov["generator"] == "lat_lon_regular"
    assert prov["source"] == "earthsci_discretizations.grids.lat_lon"


def test_dtype_f32_propagates_to_esm():
    g = lat_lon(nlon=4, nlat=4, dtype="float32")
    assert g.dtype == "float32"
    assert g.to_esm()["dtype"] == "float32"
    assert g.area.dtype == np.float32


def test_np_dtype_accepted():
    g = lat_lon(nlon=4, nlat=4, dtype=np.dtype("float32"))
    assert g.dtype == "float32"


# ------------------------------------------------------------- map_i helper --


def test_map_i_identity_when_widths_match():
    assert _map_i(0, 8, 8) == 0
    assert _map_i(7, 8, 8) == 7


def test_map_i_rounds_to_nearest():
    assert _map_i(0, 4, 8) == 1
    assert _map_i(3, 4, 8) == 7
    assert _map_i(0, 8, 4) == 0
    assert _map_i(7, 8, 4) == 3
