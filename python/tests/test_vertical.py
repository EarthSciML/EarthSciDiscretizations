"""Tests for the vertical grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §7 contract across the six coordinate kinds
(sigma, eta, z, theta, hybrid_sigma_theta, z_star): construction,
monotonicity validation, derived centers / widths, accessors, metric_eval,
and the declarative .esm round-trip shape.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

import earthsci_toolkit
from earthsci_toolkit import grids
from earthsci_toolkit.grids import VerticalGrid, vertical

# ------------------------------------------------------------------ API --


def test_exported_from_grids_namespace():
    assert grids.vertical is vertical
    assert grids.VerticalGrid is VerticalGrid


def test_coordinate_is_keyword_only():
    with pytest.raises(TypeError):
        vertical("sigma", nz=4)  # type: ignore[misc]


def test_missing_coordinate_raises():
    with pytest.raises(TypeError):
        vertical()  # type: ignore[call-arg]


def test_unknown_coordinate_rejected():
    with pytest.raises(ValueError):
        vertical(coordinate="foo", nz=4)


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", nz=4, dtype="float16")


def test_invalid_ghosts_rejected():
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", nz=4, ghosts=-1)
    with pytest.raises(TypeError):
        vertical(coordinate="sigma", nz=4, ghosts=1.5)  # type: ignore[arg-type]


def test_default_dtype_and_ghosts():
    g = vertical(coordinate="sigma", nz=4)
    assert g.dtype == "float64"
    assert g.ghosts == 0
    assert g.p0 == 1.0e5


# ------------------------------------------------------------- sigma --


def test_uniform_sigma_construction():
    g = vertical(coordinate="sigma", nz=10)
    assert isinstance(g, VerticalGrid)
    assert g.coordinate == "sigma"
    assert g.nz == 10
    assert g.levels.size == 11
    assert g.centers.size == 10
    assert g.widths.size == 10
    assert math.isclose(float(g.levels[0]), 1.0)
    assert math.isclose(float(g.levels[-1]), 0.0, abs_tol=1e-15)
    assert np.allclose(g.widths, 0.1, atol=1e-12)


def test_explicit_sigma_levels():
    lv = [1.0, 0.9, 0.7, 0.3, 0.0]
    g = vertical(coordinate="sigma", levels=lv)
    assert g.nz == 4
    assert np.allclose(g.levels, lv)
    assert np.allclose(g.centers, [0.95, 0.8, 0.5, 0.15])
    assert np.allclose(g.widths, [0.1, 0.2, 0.4, 0.3])


def test_sigma_levels_must_be_decreasing():
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", levels=[1.0, 0.5, 0.7, 0.0])


def test_sigma_levels_must_lie_in_unit_interval():
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", levels=[1.2, 0.5, 0.0])
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", levels=[1.0, 0.5, -0.1])


def test_sigma_requires_nz_or_levels():
    with pytest.raises(TypeError):
        vertical(coordinate="sigma")


def test_sigma_nz_levels_consistency():
    with pytest.raises(ValueError):
        vertical(coordinate="sigma", nz=10, levels=[1.0, 0.5, 0.0])


# --------------------------------------------------------- z / theta --


def test_z_geometric_grid():
    z = [0.0, 100.0, 300.0, 700.0, 1500.0]
    g = vertical(coordinate="z", levels=z)
    assert g.coordinate == "z"
    assert g.nz == 4
    assert np.allclose(g.levels, z)
    assert np.allclose(g.centers, [50.0, 200.0, 500.0, 1100.0])
    assert np.allclose(g.widths, [100.0, 200.0, 400.0, 800.0])


def test_z_requires_levels():
    with pytest.raises(TypeError):
        vertical(coordinate="z", nz=10)


def test_z_must_be_strictly_increasing():
    with pytest.raises(ValueError):
        vertical(coordinate="z", levels=[0.0, 100.0, 50.0])


def test_theta_grid():
    theta = [280.0, 300.0, 340.0, 400.0]
    g = vertical(coordinate="theta", levels=theta)
    assert g.coordinate == "theta"
    assert g.nz == 3
    assert np.allclose(g.centers, [290.0, 320.0, 370.0])
    assert np.allclose(g.widths, [20.0, 40.0, 60.0])


def test_z_star_grid():
    eta = [0.0, 0.25, 0.75, 1.0]
    g = vertical(coordinate="z_star", levels=eta)
    assert g.coordinate == "z_star"
    assert g.nz == 3
    assert np.allclose(g.widths, [0.25, 0.5, 0.25])


# -------------------------------------------------------------- eta --


def test_eta_hybrid_construction():
    p0 = 1.0e5
    ak = [0.0, 1000.0, 5000.0, 100.0]
    bk = [1.0, 0.7, 0.3, 0.0]
    g = vertical(coordinate="eta", ak=ak, bk=bk, p0=p0)
    assert g.coordinate == "eta"
    assert g.nz == 3
    assert g.ak.size == 4
    assert g.bk.size == 4
    # synthesized sigma = ak/p0 + bk
    expected_sigma = np.asarray(ak) / p0 + np.asarray(bk)
    assert np.allclose(g.levels, expected_sigma)


def test_eta_requires_ak_and_bk():
    with pytest.raises(TypeError):
        vertical(coordinate="eta", bk=[1.0, 0.5, 0.0])
    with pytest.raises(TypeError):
        vertical(coordinate="eta", ak=[0.0, 100.0, 200.0])


def test_eta_ak_bk_length_mismatch():
    with pytest.raises(ValueError):
        vertical(coordinate="eta", ak=[0.0, 1000.0, 100.0], bk=[1.0, 0.5, 0.25, 0.0])


def test_eta_synthesized_sigma_must_decrease():
    # bk not monotone -> synthesized sigma not monotone
    with pytest.raises(ValueError):
        vertical(
            coordinate="eta",
            ak=[0.0, 1000.0, 500.0, 100.0],
            bk=[1.0, 0.5, 0.8, 0.0],
        )


# ------------------------------------------------- hybrid_sigma_theta --


def test_hybrid_sigma_theta_uniform():
    g = vertical(coordinate="hybrid_sigma_theta", nz=5)
    assert g.coordinate == "hybrid_sigma_theta"
    assert g.nz == 5
    assert math.isclose(float(g.levels[0]), 1.0)
    assert math.isclose(float(g.levels[-1]), 0.0, abs_tol=1e-15)


def test_hybrid_sigma_theta_with_transition():
    g = vertical(coordinate="hybrid_sigma_theta", nz=4, transition=0.4)
    assert g.nz == 4


def test_hybrid_sigma_theta_transition_out_of_range():
    with pytest.raises(ValueError):
        vertical(coordinate="hybrid_sigma_theta", nz=4, transition=1.5)
    with pytest.raises(ValueError):
        vertical(coordinate="hybrid_sigma_theta", nz=4, transition=0.0)


def test_hybrid_sigma_theta_with_ak_bk():
    levels = [1.0, 0.75, 0.5, 0.25, 0.0]
    ak = [0.0, 500.0, 2000.0, 800.0, 100.0]
    bk = [1.0, 0.7, 0.4, 0.2, 0.0]
    g = vertical(
        coordinate="hybrid_sigma_theta", levels=levels, ak=ak, bk=bk, p0=1.0e5
    )
    assert g.ak.size == 5
    assert g.bk.size == 5


# -------------------------------------------------------- accessors --


def test_cell_centers_bulk_and_scalar():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0, 700.0])
    bulk = g.cell_centers()
    assert bulk.shape == (3,)
    for k, expected in enumerate([50.0, 200.0, 500.0]):
        assert math.isclose(g.cell_centers(k), expected)


def test_cell_widths_bulk_and_scalar():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0, 700.0])
    bulk = g.cell_widths()
    assert bulk.shape == (3,)
    for k, expected in enumerate([100.0, 200.0, 400.0]):
        assert math.isclose(g.cell_widths(k), expected)


def test_cell_accessors_reject_out_of_range():
    g = vertical(coordinate="sigma", nz=4)
    with pytest.raises(ValueError):
        g.cell_centers(4)
    with pytest.raises(ValueError):
        g.cell_centers(-1)
    with pytest.raises(ValueError):
        g.cell_widths(10)


def test_neighbors_interior_and_boundary():
    g = vertical(coordinate="sigma", nz=5)
    # bottom layer (k=0) has only "up"
    assert g.neighbors(0) == {"up": 1}
    # top layer (k=nz-1) has only "down"
    assert g.neighbors(4) == {"down": 3}
    # interior
    assert g.neighbors(2) == {"down": 1, "up": 3}


def test_neighbors_rejects_out_of_range():
    g = vertical(coordinate="sigma", nz=4)
    with pytest.raises(ValueError):
        g.neighbors(4)
    with pytest.raises(ValueError):
        g.neighbors(-1)


def test_n_cells_vertices_edges():
    g = vertical(coordinate="sigma", nz=7)
    assert g.n_cells == 7
    assert g.n_vertices == 8
    assert g.n_edges == 7


# ------------------------------------------------------- metric_eval --


def test_metric_eval_dz_and_z_match_accessors():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0, 700.0])
    for k in range(g.nz):
        assert math.isclose(g.metric_eval("dz", k), g.cell_widths(k))
        assert math.isclose(g.metric_eval("z", k), g.cell_centers(k))


def test_metric_eval_sigma_valid_for_sigma_like():
    g = vertical(coordinate="sigma", nz=4)
    for k in range(g.nz):
        assert math.isclose(g.metric_eval("sigma", k), g.cell_centers(k))


def test_metric_eval_sigma_rejects_z_coordinate():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0])
    with pytest.raises(ValueError):
        g.metric_eval("sigma", 0)


def test_metric_eval_pressure_for_eta():
    p0 = 1.0e5
    ak = [0.0, 1000.0, 5000.0, 100.0]
    bk = [1.0, 0.7, 0.3, 0.0]
    g = vertical(coordinate="eta", ak=ak, bk=bk, p0=p0)
    for k in range(g.nz):
        p_lo = ak[k] + bk[k] * p0
        p_hi = ak[k + 1] + bk[k + 1] * p0
        assert math.isclose(g.metric_eval("pressure", k), 0.5 * (p_lo + p_hi))


def test_metric_eval_pressure_requires_hybrid():
    g = vertical(coordinate="sigma", nz=4)
    with pytest.raises(ValueError):
        g.metric_eval("pressure", 0)


def test_metric_eval_ak_bk_averaged():
    ak = [0.0, 1000.0, 5000.0, 100.0]
    bk = [1.0, 0.7, 0.3, 0.0]
    g = vertical(coordinate="eta", ak=ak, bk=bk)
    for k in range(g.nz):
        assert math.isclose(g.metric_eval("ak", k), 0.5 * (ak[k] + ak[k + 1]))
        assert math.isclose(g.metric_eval("bk", k), 0.5 * (bk[k] + bk[k + 1]))


def test_metric_eval_rejects_unknown_name():
    g = vertical(coordinate="sigma", nz=4)
    with pytest.raises(ValueError):
        g.metric_eval("not_a_metric", 0)


def test_metric_eval_rejects_out_of_range_layer():
    g = vertical(coordinate="sigma", nz=4)
    with pytest.raises(ValueError):
        g.metric_eval("dz", 4)


# ------------------------------------------------------------ to_esm --


def test_to_esm_is_declarative_sigma():
    g = vertical(coordinate="sigma", nz=4, ghosts=2)
    doc = g.to_esm()
    assert doc["family"] == "vertical"
    assert doc["topology"] == "column"
    assert doc["ndim"] == 1
    assert doc["dtype"] == "float64"
    assert doc["ghosts"] == 2
    assert doc["n_cells"] == 4
    assert doc["n_vertices"] == 5
    assert doc["options"]["coordinate"] == "sigma"
    assert doc["options"]["nz"] == 4
    assert len(doc["options"]["levels"]) == 5
    assert "ak" not in doc["options"]
    assert "provenance" in doc


def test_to_esm_includes_hybrid_coeffs_for_eta():
    ak = [0.0, 1000.0, 5000.0, 100.0]
    bk = [1.0, 0.7, 0.3, 0.0]
    g = vertical(coordinate="eta", ak=ak, bk=bk, p0=1.0e5)
    doc = g.to_esm()
    assert doc["options"]["coordinate"] == "eta"
    assert doc["options"]["ak"] == ak
    assert doc["options"]["bk"] == bk
    assert doc["options"]["p0"] == 1.0e5


def test_to_esm_omits_derived_arrays():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0, 700.0])
    wire = json.dumps(g.to_esm())
    assert "centers" not in wire
    assert "widths" not in wire


def test_to_esm_json_roundtrip():
    g = vertical(coordinate="sigma", nz=5)
    reparsed = json.loads(json.dumps(g.to_esm()))
    assert reparsed["options"]["nz"] == 5
    assert reparsed["family"] == "vertical"


def test_provenance_carries_binding_identity():
    g = vertical(coordinate="sigma", nz=4)
    prov = g.to_esm()["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_toolkit.__version__
    assert prov["family"] == "vertical"
    assert prov["coordinate"] == "sigma"


# -------------------------------------------------------------- dtype --


def test_float32_dtype_propagates():
    g = vertical(coordinate="z", levels=[0.0, 100.0, 300.0], dtype="float32")
    assert g.dtype == "float32"
    assert g.levels.dtype == np.float32
    assert g.centers.dtype == np.float32
    assert g.widths.dtype == np.float32
    assert g.to_esm()["dtype"] == "float32"


def test_numpy_dtype_accepted():
    g = vertical(coordinate="sigma", nz=4, dtype=np.dtype("float64"))
    assert g.dtype == "float64"
