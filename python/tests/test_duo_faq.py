"""W2 (esd-un6): DUO front-door construction parity tests.

Asserts that the declarative front-door construction (``duo_subdivide_faq`` /
``primal_geometry_faq`` / ``primal_topology_faq``, all routed through the ESS
``eval_coeff`` evaluator and relational value-invention primitives) reproduces the
imperative ``duo.py`` builders:

* the integer topology (faces, vertex numbering, cell_neighbors, vertex_faces) is
  byte-identical (value-invention is integer-exact),
* the float geometry (cell_cart, lon, lat, area) agrees to libm-transcendental ULP
  (the front-door uses NumPy ``arctan2``/``arcsin``/``arccos`` while the imperative
  path uses CPython ``math.*``; both satisfy the cross-binding 1e-12 contract), and
* the edge set is identical but emitted in the ESS canonical (sorted) order rather
  than the imperative dict-insertion order.

Python mirror of the Julia ``test/test_duo_subdivision_faq.jl`` /
``test_duo_primal_geometry_faq.jl`` / topology FAQ tests.
"""

from __future__ import annotations

import importlib
import math

import numpy as np
import pytest

from earthsci_discretizations.grids.duo_primal_geometry_faq import primal_geometry_faq
from earthsci_discretizations.grids.duo_subdivide_faq import duo_subdivide_faq
from earthsci_discretizations.grids.duo_topology_faq import primal_topology_faq

# The ``grids`` package re-exports the ``duo`` *function*, shadowing the submodule;
# import the module object explicitly to reach the imperative builders.
_duo = importlib.import_module("earthsci_discretizations.grids.duo")

_LEVELS = (0, 1, 2, 3)
_R = 6.371e6
# Cross-binding contract for the DUO family (per-family override of the §4.2 1e-14
# default; L'Huilier acos/atan libm variance).
_REL_TOL = 1e-12


def _imperative_geometry(V_unit, F, R):
    """Reproduce the imperative per-cell geometry loop of ``duo()`` (pre-W2)."""
    Nc = F.shape[1]
    cell_cart = np.empty((3, Nc))
    lon = np.empty(Nc)
    lat = np.empty(Nc)
    area = np.empty(Nc)
    R2 = R * R
    for c in range(Nc):
        ai, bi, ci = int(F[0, c]), int(F[1, c]), int(F[2, c])
        a = (float(V_unit[0, ai]), float(V_unit[1, ai]), float(V_unit[2, ai]))
        b = (float(V_unit[0, bi]), float(V_unit[1, bi]), float(V_unit[2, bi]))
        cc = (float(V_unit[0, ci]), float(V_unit[1, ci]), float(V_unit[2, ci]))
        mx, my, mz = a[0] + b[0] + cc[0], a[1] + b[1] + cc[1], a[2] + b[2] + cc[2]
        n = math.sqrt(mx * mx + my * my + mz * mz)
        ux, uy, uz = mx / n, my / n, mz / n
        cell_cart[0, c], cell_cart[1, c], cell_cart[2, c] = R * ux, R * uy, R * uz
        lo, la = _duo._cart_to_lonlat(ux, uy, uz)
        lon[c], lat[c] = lo, la
        area[c] = _duo._spherical_triangle_area(a, b, cc) * R2
    return cell_cart, lon, lat, area


@pytest.mark.parametrize("level", _LEVELS)
def test_subdivide_faq_matches_imperative(level: int) -> None:
    """D3: front-door subdivision is byte-identical (vertices + faces, same order)."""
    V_faq, F_faq = duo_subdivide_faq(level)
    V_imp, F_imp = _duo._subdivide_icosahedron(level)
    assert V_faq.dtype == np.float64
    assert F_faq.dtype == np.int64
    # Encounter-order vertex numbering and face order preserved exactly.
    assert np.array_equal(F_faq, F_imp), "face table / order diverged"
    # Coordinates: same float ops (squares as products, left-assoc folds, one sqrt).
    assert np.array_equal(V_faq, V_imp), "unit-sphere vertex coords not byte-identical"


@pytest.mark.parametrize("level", _LEVELS)
def test_primal_geometry_faq_matches_imperative(level: int) -> None:
    """D2a: front-door cell geometry agrees with the imperative loop within ULP."""
    V_unit, F = _duo._subdivide_icosahedron(level)
    pg = primal_geometry_faq(V_unit, F, _R, np.float64)
    cc_i, lon_i, lat_i, area_i = _imperative_geometry(V_unit, F, _R)
    np.testing.assert_allclose(pg["cell_cart"], cc_i, rtol=_REL_TOL, atol=0.0)
    np.testing.assert_allclose(pg["lon"], lon_i, rtol=_REL_TOL, atol=1e-12)
    np.testing.assert_allclose(pg["lat"], lat_i, rtol=_REL_TOL, atol=1e-12)
    np.testing.assert_allclose(pg["area"], area_i, rtol=_REL_TOL, atol=0.0)


@pytest.mark.parametrize("level", _LEVELS)
def test_primal_topology_faq_matches_imperative(level: int) -> None:
    """D1a: value-invention topology is integer-exact; edges canonical-sorted."""
    V_unit, F = _duo._subdivide_icosahedron(level)
    Nv = V_unit.shape[1]
    topo = primal_topology_faq(F, Nv)
    edges_imp, nbr_imp = _duo._build_connectivity(F)
    vf_imp = _duo._vertex_faces(F, Nv)

    # cell_neighbors and vertex_faces are byte-identical to the imperative builders.
    assert np.array_equal(topo["cell_neighbors"], nbr_imp)
    assert topo["vertex_faces"] == vf_imp

    # Same edge SET, but emitted in canonical (sorted) order.
    faq_cols = [
        (int(topo["edges"][0, e]), int(topo["edges"][1, e]))
        for e in range(topo["edges"].shape[1])
    ]
    imp_cols = [
        (int(edges_imp[0, e]), int(edges_imp[1, e])) for e in range(edges_imp.shape[1])
    ]
    assert set(faq_cols) == set(imp_cols), "edge set diverged"
    assert faq_cols == sorted(faq_cols), "front-door edges must be canonical-sorted"
    assert len(faq_cols) == 30 * 4**level

    # edges_on_face references a valid dense edge id for every (face, slot).
    assert topo["edges_on_face"].shape == (3, F.shape[1])
    assert topo["edges_on_face"].min() >= 0
    assert topo["edges_on_face"].max() < len(faq_cols)


@pytest.mark.parametrize("level", _LEVELS)
def test_duo_grid_construction_via_front_door(level: int) -> None:
    """End-to-end: the public ``duo()`` grid matches the imperative construction."""
    g = _duo.duo(loader={"path": f"builtin://icosahedral/{level}"}, R=_R)
    V_unit, F = _duo._subdivide_icosahedron(level)
    Nv = V_unit.shape[1]
    cc_i, lon_i, lat_i, area_i = _imperative_geometry(V_unit, F, _R)
    _, nbr_imp = _duo._build_connectivity(F)
    vf_imp = _duo._vertex_faces(F, Nv)

    assert np.array_equal(g.faces, F)
    assert np.array_equal(g.vertices, _R * V_unit)
    assert np.array_equal(g.cell_neighbors, nbr_imp)
    assert g.vertex_faces == vf_imp
    np.testing.assert_allclose(g.cell_cart, cc_i, rtol=_REL_TOL, atol=0.0)
    np.testing.assert_allclose(g.lon, lon_i, rtol=_REL_TOL, atol=1e-12)
    np.testing.assert_allclose(g.lat, lat_i, rtol=_REL_TOL, atol=1e-12)
    np.testing.assert_allclose(g.area, area_i, rtol=_REL_TOL, atol=0.0)
