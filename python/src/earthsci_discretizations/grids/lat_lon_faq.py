"""Lat-lon structured-grid construction via the landed M1 elementwise FAQ path
(``eval_coeff``, :mod:`earthsci_discretizations.rules`).

Python mirror of ``src/latlon_faq.jl`` (esd-3we, S3). Declarative companion:
``discretizations/grids/latlon/rules/latlon_construction.esm``. The construction
is a semiring-FAQ (RFC ``semiring-faq-unified-ir`` §5.1/§5.2: interval index-sets
+ elementwise ``arrayop`` bodies — the M1 path, NO value-invention).

This module is the thin **evaluation bridge**: every piece of grid ARITHMETIC —
the per-row affine longitude map (``dlon = 2π/nlon``, ``lon = lon_start +
(i-½)·dlon``), the spherical-rectangle area ``R²·dlon·(sin φ_n − sin φ_s)``, and
the closed-form curvilinear metric (``g_λλ = R²cos²φ``, ``g_φφ = R²``, ``ginv``,
Jacobian ``R²|cos φ|``, and the lone non-vanishing derivative ``∂g_λλ/∂φ =
−2R²cos φ sin φ``) — flows through ESS's single AST pathway (``eval_coeff``). The
structural integer arrays (ragged-row-major flattening, periodic-longitude
neighbor linearization, the reduced-Gaussian nearest-center N/S remap, boundary
masks) are pure index logic. Byte/ULP-identical to the imperative ``LatLonGrid``
accessors and to the committed Julia-reference
``tests/conformance/grids/latlon/construction/golden.json``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..rules import eval_coeff

__all__ = ["LatLonConstruction", "lat_lon_construction_faq"]


def _mk(op: str, *args) -> dict:
    return {"op": op, "args": list(args)}


# π inlined as the binary64 literal 2π (identical to the imperative ``2 * pi``).
_DLON = _mk("/", 6.283185307179586, "n")  # 2π / n
_LON_CENTER = _mk("+", "lon_start", _mk("*", _mk("-", "i", 0.5), "dlon"))
_LAT_WIDTH = _mk("-", "e_hi", "e_lo")
# area = ((R*R)*dlon) * (sin(lat_n) - sin(lat_s))
_AREA = _mk(
    "*",
    _mk("*", _mk("*", "R", "R"), "dlon"),
    _mk("-", _mk("sin", "lat_n"), _mk("sin", "lat_s")),
)
_COS = _mk("cos", "lat")
_GLL = _mk("*", _mk("*", _mk("*", "R", "R"), "cos_lat"), "cos_lat")  # R²cos²φ
_GPP = _mk("*", "R", "R")  # R²
_INV = _mk("/", 1.0, "x")
_JAC = _mk("*", _mk("*", "R", "R"), _mk("abs", _mk("cos", "lat")))  # R²|cosφ|
# ∂g_λλ/∂φ = -2*R*R*cos(φ)*sin(φ)
_DGLL = _mk(
    "*",
    _mk("*", _mk("*", _mk("*", -2.0, "R"), "R"), _mk("cos", "phi")),
    _mk("sin", "phi"),
)


@dataclass(frozen=True)
class LatLonConstruction:
    """Materialized lat-lon construction (mirrors the imperative accessors and
    the ``construction/golden.json`` contract). Per-cell arrays use the flat
    ragged-row-major layout (cell ``(j, i)`` at flat index ``row_off[j] + i``).
    The metric is stored as the golden's named scalar components."""

    cell_centers_lon: np.ndarray
    cell_centers_lat: np.ndarray
    cell_widths_lon: np.ndarray
    cell_widths_lat: np.ndarray
    cell_volume: np.ndarray
    lat_edges: np.ndarray
    lat_centers: np.ndarray
    metric_g_lonlon: np.ndarray
    metric_g_latlat: np.ndarray
    metric_ginv_lonlon: np.ndarray
    metric_jacobian: np.ndarray
    dg_lonlon_dlat: np.ndarray
    neighbor_lon_minus: np.ndarray
    neighbor_lon_plus: np.ndarray
    neighbor_lat_minus: np.ndarray
    neighbor_lat_plus: np.ndarray
    boundary_lat_lower: np.ndarray
    boundary_lat_upper: np.ndarray


def _latlon_map_i(i: int, from_n: int, to_n: int) -> int:
    """Reduced-Gaussian nearest-center N/S column remap (0-based), mirroring the
    imperative ``_latlon_map_i`` (the relational rank/join materialized
    structurally)."""
    if from_n == to_n:
        return i
    frac = (i + 0.5) / from_n
    k = math.floor(frac * to_n)
    return min(max(k, 0), to_n - 1)


def lat_lon_construction_faq(grid) -> LatLonConstruction:
    """Materialize a ``LatLonGrid``'s construction from the declarative FAQ,
    routing all arithmetic through the landed ESS AST evaluator (``eval_coeff``).
    """
    np_dtype = np.dtype(grid.dtype).type
    nlat = int(grid.nlat)
    per_row = [int(x) for x in grid.nlon_per_row]
    nc = int(sum(per_row))
    r = float(grid.R)
    lon_start = float(grid.lon_start)
    lat_edges = np.asarray(grid.lat_edges, dtype=np_dtype)
    lat_centers = np.asarray(grid.lat_centers, dtype=np_dtype)

    row_off = [0] * (nlat + 1)
    for j in range(nlat):
        row_off[j + 1] = row_off[j] + per_row[j]

    cc_lon = np.empty(nc, dtype=np_dtype)
    cc_lat = np.empty(nc, dtype=np_dtype)
    cw_lon = np.empty(nc, dtype=np_dtype)
    cw_lat = np.empty(nc, dtype=np_dtype)
    cell_volume = np.empty(nc, dtype=np_dtype)
    g_lonlon = np.empty(nc, dtype=np_dtype)
    g_latlat = np.empty(nc, dtype=np_dtype)
    ginv_lonlon = np.empty(nc, dtype=np_dtype)
    jacobian = np.empty(nc, dtype=np_dtype)
    dg = np.empty(nc, dtype=np_dtype)

    # --- per-row elementwise FAQ (affine lon, trig metric, area) ------------
    for j in range(nlat):
        n_i = per_row[j]
        base = row_off[j]
        dlon = eval_coeff(_DLON, {"n": float(n_i)})
        lat_c = float(lat_centers[j])
        lat_s = float(lat_edges[j])
        lat_n = float(lat_edges[j + 1])
        lat_w = eval_coeff(_LAT_WIDTH, {"e_hi": lat_n, "e_lo": lat_s})
        area = eval_coeff(_AREA, {"R": r, "dlon": dlon, "lat_n": lat_n, "lat_s": lat_s})
        cos_lat = eval_coeff(_COS, {"lat": lat_c})
        g_ll = eval_coeff(_GLL, {"R": r, "cos_lat": cos_lat})
        g_pp = eval_coeff(_GPP, {"R": r})
        inv_ll = eval_coeff(_INV, {"x": g_ll}) if g_ll > 0.0 else math.inf
        jac = eval_coeff(_JAC, {"R": r, "lat": lat_c})
        dgll = eval_coeff(_DGLL, {"R": r, "phi": lat_c})
        for i in range(n_i):
            k = base + i
            cc_lon[k] = eval_coeff(
                _LON_CENTER, {"lon_start": lon_start, "i": float(i + 1), "dlon": dlon}
            )
            cc_lat[k] = lat_c
            cw_lon[k] = dlon
            cw_lat[k] = lat_w
            cell_volume[k] = area
            g_lonlon[k] = g_ll
            g_latlat[k] = g_pp
            ginv_lonlon[k] = inv_ll
            jacobian[k] = jac
            dg[k] = dgll

    # --- neighbor maps (flat ragged-row-major linear ids) -------------------
    nbr_lon_minus = np.empty(nc, dtype=np.int64)
    nbr_lon_plus = np.empty(nc, dtype=np.int64)
    for j in range(nlat):
        n_i = per_row[j]
        base = row_off[j]
        for i in range(n_i):
            nbr_lon_minus[base + i] = base + (i + n_i - 1) % n_i
            nbr_lon_plus[base + i] = base + (i + 1) % n_i

    def _nbr_lat(offset: int) -> np.ndarray:
        out = np.empty(nc, dtype=np.int64)
        for j in range(nlat):
            n_i = per_row[j]
            base = row_off[j]
            jj = j + offset
            if 0 <= jj < nlat:
                n_t = per_row[jj]
                base_t = row_off[jj]
                for i in range(n_i):
                    out[base + i] = base_t + _latlon_map_i(i, n_i, n_t)
            else:
                for i in range(n_i):
                    out[base + i] = -1  # pole: pole_policy = none
        return out

    neighbor_lat_minus = _nbr_lat(-1)
    neighbor_lat_plus = _nbr_lat(+1)

    # --- boundary masks (lon wraps → all-false; lat marks first/last row) ---
    def _bd_lat(target: int) -> np.ndarray:
        out = np.zeros(nc, dtype=bool)
        for j in range(nlat):
            n_i = per_row[j]
            base = row_off[j]
            if j == target:
                out[base : base + n_i] = True
        return out

    boundary_lat_lower = _bd_lat(0)
    boundary_lat_upper = _bd_lat(nlat - 1)

    return LatLonConstruction(
        cell_centers_lon=cc_lon,
        cell_centers_lat=cc_lat,
        cell_widths_lon=cw_lon,
        cell_widths_lat=cw_lat,
        cell_volume=cell_volume,
        lat_edges=lat_edges.copy(),
        lat_centers=lat_centers.copy(),
        metric_g_lonlon=g_lonlon,
        metric_g_latlat=g_latlat,
        metric_ginv_lonlon=ginv_lonlon,
        metric_jacobian=jacobian,
        dg_lonlon_dlat=dg,
        neighbor_lon_minus=nbr_lon_minus,
        neighbor_lon_plus=nbr_lon_plus,
        neighbor_lat_minus=neighbor_lat_minus,
        neighbor_lat_plus=neighbor_lat_plus,
        boundary_lat_lower=boundary_lat_lower,
        boundary_lat_upper=boundary_lat_upper,
    )
