"""Tests for the cartesian grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §7 contract: 1D/2D/3D uniform + non-uniform
construction, accessor semantics (cell_centers / cell_widths / cell_volume /
neighbors / metric_eval), determinism, error contract (§9), and declarative
.esm lowering shape (including byte-diff against committed fixtures).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import earthsci_discretizations
from earthsci_discretizations import grids
from earthsci_discretizations.grids import CartesianGrid, cartesian

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "discretizations" / "grids" / "cartesian"


# ---------------------------------------------------------------- exports --


def test_exported_from_grids_namespace():
    assert grids.cartesian is cartesian
    assert grids.CartesianGrid is CartesianGrid


# ---------------------------------------------------- construction: uniform


def test_uniform_1d_construction():
    g = cartesian(nx=10, extent=[(0.0, 1.0)])
    assert isinstance(g, CartesianGrid)
    assert g.ndim == 1
    assert g.n == (10,)
    assert g.ghosts == 0
    assert g.uniform == (True,)
    assert g.family == "cartesian"
    assert g.topology == "rectilinear"
    assert g.dtype == "float64"
    assert g.edges[0].shape == (11,)
    assert g.centers[0].shape == (10,)
    assert float(g.edges[0][0]) == pytest.approx(0.0)
    assert float(g.edges[0][-1]) == pytest.approx(1.0)
    np.testing.assert_allclose(g.widths[0], 0.1, atol=1e-12)


def test_uniform_2d_vector_of_pairs_extent():
    g = cartesian(nx=4, ny=8, extent=[(0.0, 2.0), (0.0, 4.0)])
    assert g.ndim == 2
    assert g.n == (4, 8)
    assert g.extent == ((0.0, 2.0), (0.0, 4.0))
    np.testing.assert_allclose(g.widths[0], 0.5, atol=1e-12)
    np.testing.assert_allclose(g.widths[1], 0.5, atol=1e-12)


def test_uniform_2d_matrix_extent_ndarray():
    ext = np.array([[0.0, 0.0], [2.0, 4.0]])  # shape (2, 2)
    g = cartesian(nx=4, ny=8, extent=ext)
    assert g.extent == ((0.0, 2.0), (0.0, 4.0))


def test_uniform_3d_construction():
    g = cartesian(nx=2, ny=3, nz=4, extent=[(0.0, 1.0), (0.0, 1.5), (0.0, 2.0)])
    assert g.ndim == 3
    assert g.n == (2, 3, 4)
    assert g.uniform == (True, True, True)
    assert g.n_cells == 24


def test_uniform_matrix_extent_2x3_for_3d():
    ext = np.array([[0.0, 0.0, 0.0], [1.0, 1.5, 2.0]])  # shape (2, 3)
    g = cartesian(nx=2, ny=3, nz=4, extent=ext)
    assert g.extent == ((0.0, 1.0), (0.0, 1.5), (0.0, 2.0))


def test_dtype_float32():
    g = cartesian(nx=4, extent=[(0.0, 1.0)], dtype="float32")
    assert g.dtype == "float32"
    assert g.centers[0].dtype == np.float32
    assert g.edges[0].dtype == np.float32


def test_dtype_numpy_dtype_accepted():
    g = cartesian(nx=4, extent=[(0.0, 1.0)], dtype=np.dtype("float32"))
    assert g.dtype == "float32"


def test_ghosts_option():
    g = cartesian(nx=4, extent=[(0.0, 1.0)], ghosts=3)
    assert g.ghosts == 3


# ------------------------------------------- construction: non-uniform --


def test_nonuniform_2d_edges():
    xe = [0.0, 0.1, 0.3, 0.7, 1.0]
    ye = [-1.0, 0.0, 0.5, 1.0]
    g = cartesian(edges=[xe, ye])
    assert g.ndim == 2
    assert g.n == (4, 3)
    assert g.uniform == (False, False)
    assert g.extent == ((0.0, 1.0), (-1.0, 1.0))
    np.testing.assert_allclose(g.widths[0], [0.1, 0.2, 0.4, 0.3])
    np.testing.assert_allclose(g.widths[1], [1.0, 0.5, 0.5])


def test_nonuniform_detects_uniform_widths():
    # Edges look non-uniform but widths are uniform → uniform flag True.
    g = cartesian(edges=[[0.0, 0.25, 0.5, 0.75, 1.0]])
    assert g.uniform == (True,)


# ----------------------------------------------------------- accessors --


def test_cell_centers_1d_uniform():
    g = cartesian(nx=4, extent=[(0.0, 1.0)])
    # Centers at 0.125, 0.375, 0.625, 0.875
    assert g.cell_centers(0)[0] == pytest.approx(0.125)
    assert g.cell_centers(1)[0] == pytest.approx(0.375)
    assert g.cell_centers(2)[0] == pytest.approx(0.625)
    assert g.cell_centers(3)[0] == pytest.approx(0.875)


def test_cell_centers_2d():
    g = cartesian(nx=2, ny=2, extent=[(0.0, 2.0), (10.0, 14.0)])
    assert g.cell_centers(0, 0) == (0.5, 11.0)
    assert g.cell_centers(1, 1) == (1.5, 13.0)


def test_cell_widths_2d_uniform():
    g = cartesian(nx=4, ny=2, extent=[(0.0, 1.0), (0.0, 0.5)])
    assert g.cell_widths(0, 0) == pytest.approx((0.25, 0.25))
    assert g.cell_widths(3, 1) == pytest.approx((0.25, 0.25))


def test_cell_volume_3d():
    g3 = cartesian(nx=2, ny=2, nz=2, extent=[(0.0, 2.0), (0.0, 4.0), (0.0, 8.0)])
    assert g3.cell_volume(0, 0, 0) == pytest.approx(1.0 * 2.0 * 4.0)


def test_cell_volume_nonuniform_1d():
    g_nu = cartesian(edges=[[0.0, 0.5, 1.5]])  # widths 0.5, 1.0
    assert g_nu.cell_volume(0) == pytest.approx(0.5)
    assert g_nu.cell_volume(1) == pytest.approx(1.0)


def test_neighbors_interior_and_corner_cells():
    g = cartesian(nx=3, ny=3, extent=[(0.0, 1.0), (0.0, 1.0)])
    # Interior cell (1, 1) has 4 neighbors
    nb = g.neighbors(1, 1)
    assert len(nb) == 4
    assert nb[(0, -1)] == (0, 1)
    assert nb[(0, +1)] == (2, 1)
    assert nb[(1, -1)] == (1, 0)
    assert nb[(1, +1)] == (1, 2)
    # Corner cell (0, 0) has 2 neighbors
    nb_corner = g.neighbors(0, 0)
    assert len(nb_corner) == 2
    assert (0, +1) in nb_corner
    assert (1, +1) in nb_corner
    assert (0, -1) not in nb_corner
    assert (1, -1) not in nb_corner


def test_neighbors_1d_boundaries():
    g = cartesian(nx=4, extent=[(0.0, 1.0)])
    left = g.neighbors(0)
    right = g.neighbors(3)
    assert list(left.keys()) == [(0, +1)]
    assert list(right.keys()) == [(0, -1)]
    # Middle cell
    mid = g.neighbors(1)
    assert len(mid) == 2


def test_metric_eval_volume_dx_dy_face_areas_2d():
    g = cartesian(nx=2, ny=4, extent=[(0.0, 2.0), (0.0, 4.0)])
    # dx = 1.0, dy = 1.0 → volume = 1.0
    assert g.metric_eval("volume", 0, 0) == pytest.approx(1.0)
    assert g.metric_eval("dx", 0, 0) == pytest.approx(1.0)
    assert g.metric_eval("dy", 0, 0) == pytest.approx(1.0)
    # 2D face_area_x = dy, face_area_y = dx
    assert g.metric_eval("face_area_x", 0, 0) == pytest.approx(1.0)
    assert g.metric_eval("face_area_y", 0, 0) == pytest.approx(1.0)
    assert g.metric_eval("jacobian", 0, 0) == 1.0
    gI = g.metric_eval("g", 0, 0)
    assert gI == ((1.0, 0.0), (0.0, 1.0))


def test_metric_eval_face_area_3d():
    g = cartesian(nx=2, ny=2, nz=2, extent=[(0.0, 2.0), (0.0, 4.0), (0.0, 8.0)])
    # widths: 1.0, 2.0, 4.0 → face_area_x = dy*dz, face_area_y = dx*dz, face_area_z = dx*dy
    assert g.metric_eval("face_area_x", 0, 0, 0) == pytest.approx(2.0 * 4.0)
    assert g.metric_eval("face_area_y", 0, 0, 0) == pytest.approx(1.0 * 4.0)
    assert g.metric_eval("face_area_z", 0, 0, 0) == pytest.approx(1.0 * 2.0)
    gI = g.metric_eval("g", 0, 0, 0)
    assert gI == (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


def test_metric_eval_rejects_unknown_axis_or_name():
    g1 = cartesian(nx=2, extent=[(0.0, 1.0)])
    with pytest.raises(ValueError):
        g1.metric_eval("dy", 0)
    with pytest.raises(ValueError):
        g1.metric_eval("dz", 0)
    with pytest.raises(ValueError):
        g1.metric_eval("bogus", 0)


# -------------------------------------------------- topology consistency --


@pytest.mark.parametrize("nx,ext", [(1, (0.0, 1.0)), (5, (0.0, 1.0)), (16, (-3.0, 7.5))])
def test_topology_edges_bound_centers_1d(nx, ext):
    g = cartesian(nx=nx, extent=[ext])
    for i in range(nx):
        assert float(g.edges[0][i]) < float(g.centers[0][i]) < float(g.edges[0][i + 1])
    assert float(np.sum(g.widths[0])) == pytest.approx(ext[1] - ext[0])
    assert float(g.edges[0][0]) == pytest.approx(ext[0])
    assert float(g.edges[0][-1]) == pytest.approx(ext[1])


def test_topology_nonuniform_partition_of_unity_2d():
    g = cartesian(edges=[[0.0, 0.3, 1.0], [0.0, 0.5, 2.0]])
    total = 0.0
    for i in range(g.n[0]):
        for j in range(g.n[1]):
            total += g.cell_volume(i, j)
    assert total == pytest.approx(1.0 * 2.0)


# --------------------------------------------------------- metric accuracy --


def test_metric_accuracy_uniform_dx_exact():
    nx = 7
    ext = (-1.0, 1.0)
    g = cartesian(nx=nx, extent=[ext])
    expected_dx = (ext[1] - ext[0]) / nx
    for i in range(nx):
        assert g.metric_eval("dx", i) == pytest.approx(expected_dx)


def test_metric_accuracy_nonuniform_widths_exact():
    xe = [0.0, 1.0, 1.5, 4.0, 4.25]
    g = cartesian(edges=[xe])
    for i in range(len(xe) - 1):
        assert g.metric_eval("dx", i) == pytest.approx(xe[i + 1] - xe[i])


# ------------------------------------------------- error contract (§9) --


def test_missing_required_options_raise_typeerror():
    with pytest.raises(TypeError):
        cartesian()
    with pytest.raises(TypeError):
        cartesian(nx=4)  # missing extent
    with pytest.raises(TypeError):
        cartesian(nx=4, ny=2)  # 2D missing extent
    with pytest.raises(TypeError):
        # 3D with nx and nz but no ny
        cartesian(nx=4, nz=2, extent=[(0.0, 1.0), (0.0, 1.0)])


def test_keyword_only_api():
    with pytest.raises(TypeError):
        cartesian(4, [(0.0, 1.0)])  # type: ignore[misc]


def test_invalid_extent_or_negative_ghosts_raise():
    with pytest.raises(ValueError):
        cartesian(nx=4, extent=[(1.0, 0.0)])  # hi <= lo
    with pytest.raises(ValueError):
        cartesian(nx=0, extent=[(0.0, 1.0)])  # n < 1
    with pytest.raises(ValueError):
        cartesian(nx=4, extent=[(0.0, 1.0)], ghosts=-1)


def test_extent_shape_mismatch_raises():
    # 3D grid needs length-3 extent, not length-2
    with pytest.raises(ValueError):
        cartesian(nx=2, ny=2, nz=2, extent=[(0.0, 1.0), (0.0, 1.0)])
    # 2D matrix form with wrong column count
    with pytest.raises(ValueError):
        cartesian(nx=2, ny=2, extent=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]))


def test_nonuniform_edges_strictly_increasing():
    with pytest.raises(ValueError):
        cartesian(edges=[[0.0, 0.5, 0.5, 1.0]])  # duplicate
    with pytest.raises(ValueError):
        cartesian(edges=[[0.0, 0.5, 0.25, 1.0]])  # decreasing
    with pytest.raises(ValueError):
        cartesian(edges=[[0.0]])  # too short


def test_cannot_mix_edges_and_nx_or_extent():
    with pytest.raises(ValueError):
        cartesian(nx=4, edges=[[0.0, 1.0, 2.0]])
    with pytest.raises(ValueError):
        cartesian(extent=[(0.0, 1.0)], edges=[[0.0, 1.0, 2.0]])


def test_nonuniform_ndim_cap_le_3():
    with pytest.raises(ValueError):
        cartesian(
            edges=[[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
        )


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        cartesian(nx=4, extent=[(0.0, 1.0)], dtype="float16")


def test_bool_not_accepted_for_int_fields():
    with pytest.raises(TypeError):
        cartesian(nx=True, extent=[(0.0, 1.0)])  # type: ignore[arg-type]


def test_index_bounds_checked():
    g = cartesian(nx=3, extent=[(0.0, 1.0)])
    with pytest.raises(ValueError):
        g.cell_centers(3)
    with pytest.raises(ValueError):
        g.cell_centers(-1)
    with pytest.raises(TypeError):
        g.cell_centers(0, 0)  # too many indices
    with pytest.raises(TypeError):
        g.cell_centers()


# ------------------------------------------------------------- to_esm --


def test_to_esm_uniform_fields():
    g = cartesian(nx=4, ny=2, extent=[(0.0, 1.0), (0.0, 0.5)], ghosts=2)
    esm = g.to_esm()
    assert esm["family"] == "cartesian"
    assert esm["topology"] == "rectilinear"
    assert esm["dtype"] == "float64"
    assert esm["ndim"] == 2
    assert esm["ghosts"] == 2
    assert esm["n_cells"] == 4 * 2
    assert esm["n"] == [4, 2]
    assert esm["uniform"] == [True, True]
    assert esm["extent"][0] == [0.0, 1.0]
    assert esm["extent"][1] == [0.0, 0.5]
    # Uniform: no explicit edges
    assert "edges" not in esm
    prov = esm["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_discretizations.__version__
    assert prov["source"] == "earthsci_discretizations.grids.cartesian"


def test_to_esm_nonuniform_includes_edges():
    xe = [0.0, 0.1, 0.4, 1.0]
    g = cartesian(edges=[xe])
    esm = g.to_esm()
    assert esm["uniform"] == [False]
    assert esm["edges"][0] == xe


def test_to_esm_dtype_reflects_element_type():
    g32 = cartesian(nx=2, extent=[(0.0, 1.0)], dtype="float32")
    assert g32.to_esm()["dtype"] == "float32"


def test_to_esm_is_json_serializable():
    g = cartesian(nx=4, ny=2, extent=[(0.0, 1.0), (0.0, 0.5)])
    s = json.dumps(g.to_esm())
    assert isinstance(s, str)
    assert '"family": "cartesian"' in s


# ------------------------- determinism: same opts → equal .esm payloads --


def test_determinism_to_esm_equal_for_same_opts():
    opts = dict(nx=8, ny=8, extent=[(0.0, 1.0), (0.0, 1.0)], ghosts=1)
    a = cartesian(**opts).to_esm()
    b = cartesian(**opts).to_esm()
    assert a == b


# ---------------- cross-binding fixture round-trip (byte-diff, no prov) --


def _strip_provenance(esm: dict) -> dict:
    d = dict(esm)
    d.pop("provenance", None)
    return d


def _fixture_cases():
    # (filename, opts) pairs mirroring discretizations/grids/cartesian/README.md
    return [
        ("uniform_1d_n16.esm", dict(nx=16, extent=[(0.0, 1.0)])),
        ("uniform_1d_n64.esm", dict(nx=64, extent=[(0.0, 1.0)])),
        ("uniform_1d_n256.esm", dict(nx=256, extent=[(0.0, 1.0)])),
        ("uniform_2d_n64.esm", dict(nx=64, ny=64, extent=[(0.0, 1.0), (0.0, 1.0)])),
        (
            "uniform_3d_n16.esm",
            dict(nx=16, ny=16, nz=16, extent=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]),
        ),
        (
            "nonuniform_2d_stretched.esm",
            dict(edges=[[0.0, 0.1, 0.3, 0.7, 1.0], [-1.0, 0.0, 0.5, 1.0]]),
        ),
    ]


@pytest.mark.parametrize("filename,opts", _fixture_cases())
def test_fixture_round_trip_minus_provenance(filename, opts):
    fixture_path = FIXTURE_DIR / filename
    if not fixture_path.exists():
        pytest.skip(f"fixture {filename} not present")
    expected = json.loads(fixture_path.read_text())
    g = cartesian(**opts)
    actual = _strip_provenance(g.to_esm())
    # Fixture is already provenance-stripped per README.
    assert actual == expected


# ------------------------------------------------------------ repr --


def test_repr_round_trips_key_fields():
    g = cartesian(nx=4, extent=[(0.0, 1.0)], ghosts=2)
    s = repr(g)
    assert "CartesianGrid" in s
    assert "ndim=1" in s
    assert "ghosts=2" in s
