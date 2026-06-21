"""DUO icosahedral primal-cell geometry via the landed ESS AST evaluator
(``eval_coeff``, :mod:`earthsci_discretizations.rules`).

Python mirror of ``src/primal_geometry_faq.jl`` (D2a / esd-heg.6). Declarative
companion: ``discretizations/grids/duo/faq/primal_geometry.esm``. That document
expresses the per-cell geometry of the DUO builder — the cell-center centroid
(normalized mean of the three corner unit vectors), the geographic ``lon``/``lat``
(``atan2``/``asin``), and the spherical-excess cell area (L'Huilier) — as the RFC
``semiring-faq-unified-ir`` §8.1 geometry-FAQ pattern.

This file is the thin **evaluation bridge**: it shapes the unit vertices + faces
into the per-cell scalar inputs and routes every piece of grid ARITHMETIC through
ESS's single AST pathway (``eval_coeff``). ESD hosts no shadow evaluator
(AGENTS.md single-pathway invariant); the determinism / cross-binding byte-identity
of the output lives entirely in ``earthsci_toolkit``. The only ESD-side logic is
gathering the corner coordinates and the ``float`` coercions.

The transcription mirrors the imperative per-cell loop (:mod:`duo`) bit-for-bit
(``float64``): squares are products (``x*x``, not ``**``); the centroid sum, the
squared norm, and the three side-arc dot products fold in corner / space order;
the L'Huilier tan-product is left-folded; ``area = (4·atan(√t)) · (R·R)`` keeps the
imperative association ``a_unit * R2``. The single intentional divergence is
dropping the ``clamp`` / ``max`` clipping guards: on a valid icosahedral mesh every
``acos`` / ``asin`` argument and the L'Huilier radicand already lie in range, so the
guards never fire and the bytes match.
"""

from __future__ import annotations

import numpy as np

from ..rules import eval_coeff

__all__ = ["primal_geometry_faq"]


# AST nodes mirroring the equation bodies of
# ``discretizations/grids/duo/faq/primal_geometry.esm``. Held as module consts so
# per-cell evaluation does not rebuild the trees. Scalar corner leaves: a1..a3
# (corner 1 = face vertex 1), b1..b3 (corner 2), c1..c3 (corner 3), and radius R.
def _pg_mk(op: str, *args) -> dict:
    return {"op": op, "args": list(args)}


# centroid: m = a+b+c per component; n = sqrt(mx*mx+my*my+mz*mz); u = m/n.
_PG_MX = _pg_mk("+", "a1", "b1", "c1")
_PG_MY = _pg_mk("+", "a2", "b2", "c2")
_PG_MZ = _pg_mk("+", "a3", "b3", "c3")
_PG_NSQ = _pg_mk(
    "+",
    _pg_mk("*", _PG_MX, _PG_MX),
    _pg_mk("*", _PG_MY, _PG_MY),
    _pg_mk("*", _PG_MZ, _PG_MZ),
)  # squares as *, not **
_PG_NRM = _pg_mk("sqrt", _PG_NSQ)
_PG_UX = _pg_mk("/", _PG_MX, _PG_NRM)
_PG_UY = _pg_mk("/", _PG_MY, _PG_NRM)
_PG_UZ = _pg_mk("/", _PG_MZ, _PG_NRM)

# cell-center cartesian (R-scaled), geographic lon/lat of the unit centroid.
_PG_CART_X = _pg_mk("*", "R", _PG_UX)
_PG_CART_Y = _pg_mk("*", "R", _PG_UY)
_PG_CART_Z = _pg_mk("*", "R", _PG_UZ)
_PG_LON = _pg_mk("atan2", _PG_UY, _PG_UX)
_PG_LAT = _pg_mk("asin", _PG_UZ)  # no clamp

# area (L'Huilier): da=acos(b·c), db=acos(c·a), dc=acos(a·b); s=(da+db+dc)/2;
# t=tan(s/2)·tan((s-da)/2)·tan((s-db)/2)·tan((s-dc)/2); area = (4·atan(√t)) · (R·R).
_PG_DOT_BC = _pg_mk("+", _pg_mk("*", "b1", "c1"), _pg_mk("*", "b2", "c2"), _pg_mk("*", "b3", "c3"))
_PG_DOT_CA = _pg_mk("+", _pg_mk("*", "c1", "a1"), _pg_mk("*", "c2", "a2"), _pg_mk("*", "c3", "a3"))
_PG_DOT_AB = _pg_mk("+", _pg_mk("*", "a1", "b1"), _pg_mk("*", "a2", "b2"), _pg_mk("*", "a3", "b3"))
_PG_DA = _pg_mk("acos", _PG_DOT_BC)
_PG_DB = _pg_mk("acos", _PG_DOT_CA)
_PG_DC = _pg_mk("acos", _PG_DOT_AB)
_PG_S = _pg_mk("*", 0.5, _pg_mk("+", _PG_DA, _PG_DB, _PG_DC))


def _pg_half(x: dict) -> dict:
    return _pg_mk("/", x, 2.0)


_PG_T = _pg_mk(
    "*",
    _pg_mk("tan", _pg_half(_PG_S)),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DA))),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DB))),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DC))),
)
_PG_AREA = _pg_mk(
    "*",
    _pg_mk("*", 4.0, _pg_mk("atan", _pg_mk("sqrt", _PG_T))),
    _pg_mk("*", "R", "R"),
)  # no max guard


def primal_geometry_faq(
    vertices_unit: np.ndarray,
    faces: np.ndarray,
    R: float,
    np_dtype: type,
) -> dict[str, np.ndarray]:
    """Materialize the DUO primal-cell geometry from the declarative FAQ, routing
    all arithmetic through the landed ESS AST evaluator (``eval_coeff``).

    Parameters
    ----------
    vertices_unit : np.ndarray
        ``(3, Nv)`` cartesian **unit-sphere** vertex coordinates (NOT R-scaled),
        as produced by :func:`duo_subdivide_faq`.
    faces : np.ndarray
        ``(3, Nc)`` 0-based vertex indices per triangular cell.
    R : float
        Sphere radius.
    np_dtype : type
        Output element dtype (``np.float64`` or ``np.float32``).

    Returns
    -------
    dict with keys ``cell_cart`` ``(3, Nc)`` R-scaled cell-center cartesian,
    ``lon`` ``(Nc,)`` longitude (rad), ``lat`` ``(Nc,)`` latitude (rad), and
    ``area`` ``(Nc,)`` spherical-triangle area (R² scaled).

    Byte-identical (``float64``) to the imperative cell-geometry loop (:mod:`duo`):
    the per-cell bindings and ASTs mirror the imperative float ops exactly, and the
    ``clamp`` / ``max`` guards are dropped because they never fire on a valid mesh.
    """
    V = vertices_unit
    F = faces
    Nc = int(F.shape[1])
    Rf = float(R)

    cell_cart = np.empty((3, Nc), dtype=np_dtype)
    lon = np.empty(Nc, dtype=np_dtype)
    lat = np.empty(Nc, dtype=np_dtype)
    area = np.empty(Nc, dtype=np_dtype)

    for c in range(Nc):
        a_i = int(F[0, c])
        b_i = int(F[1, c])
        c_i = int(F[2, c])
        # One bindings dict per cell drives every co-derived metric (centroid →
        # cart/lon/lat, and the L'Huilier area), mirroring the single imperative loop.
        b = {
            "a1": float(V[0, a_i]), "a2": float(V[1, a_i]), "a3": float(V[2, a_i]),
            "b1": float(V[0, b_i]), "b2": float(V[1, b_i]), "b3": float(V[2, b_i]),
            "c1": float(V[0, c_i]), "c2": float(V[1, c_i]), "c3": float(V[2, c_i]),
            "R": Rf,
        }
        cell_cart[0, c] = eval_coeff(_PG_CART_X, b)
        cell_cart[1, c] = eval_coeff(_PG_CART_Y, b)
        cell_cart[2, c] = eval_coeff(_PG_CART_Z, b)
        lon[c] = eval_coeff(_PG_LON, b)
        lat[c] = eval_coeff(_PG_LAT, b)
        area[c] = eval_coeff(_PG_AREA, b)

    return {"cell_cart": cell_cart, "lon": lon, "lat": lat, "area": area}
