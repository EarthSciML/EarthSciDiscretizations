"""Vertical (1-D column) structured-grid construction via the landed M1
elementwise FAQ path (``eval_coeff``, :mod:`earthsci_discretizations.rules`).

Python mirror of ``src/vertical_faq.jl`` (esd-3we, S2). Declarative companion:
``discretizations/grids/vertical/rules/vertical_construction.esm``. The
construction is a semiring-FAQ (RFC ``semiring-faq-unified-ir`` §5.1/§5.2: two
interval index-sets — ``interfaces`` (nz+1 half-levels) and ``layers`` (nz cells)
— plus elementwise ``arrayop`` bodies; the M1 path, NO value-invention).

This module is the thin **evaluation bridge**: every piece of vertical grid
ARITHMETIC — the per-coordinate level synthesis, the cell-center/width
derivations, the cell-volume, and the named layer metrics — flows through ESS's
single AST pathway (``eval_coeff``). The structural integer arrays (±k neighbor
maps with the ``-1`` off-column sentinel, boundary masks) are pure index logic.
Byte-identical to the imperative ``VerticalGrid`` accessors and to the committed
Julia-reference ``tests/conformance/grids/vertical/construction/golden.json``
(strict 0.0 tolerance — no transcendentals).

Level synthesis mirrors the imperative generator branch-for-branch: ``eta`` is
the divide-add ``ak[k]/p0 + bk[k]``; uniform ``sigma`` / ``hybrid_sigma_theta``
is the affine ``1 - k/nz`` (used only when it reproduces the stored levels
bit-for-bit; an explicit-level grid carries them as DATA); ``z`` / ``theta`` /
``z_star`` supply their levels verbatim.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..rules import eval_coeff

__all__ = ["VerticalConstruction", "vertical_construction_faq"]


def _mk(op: str, *args) -> dict:
    return {"op": op, "args": list(args)}


_SIGMA = _mk("-", 1, _mk("/", "k", "nz"))  # 1 - k/nz  (k = 0-based interface)
_ETA = _mk("+", _mk("/", "ak", "p0"), "bk")  # ak/p0 + bk
_AVG2 = _mk("/", _mk("+", "a", "b"), 2)  # (a + b) / 2
_WIDTH = _mk("abs", _mk("-", "b", "a"))  # abs(b - a)
_PRESSURE = _mk(
    "/",
    _mk(
        "+",
        _mk("+", "ak_lo", _mk("*", "bk_lo", "p0")),
        _mk("+", "ak_hi", _mk("*", "bk_hi", "p0")),
    ),
    2,
)  # ((ak_lo + bk_lo*p0) + (ak_hi + bk_hi*p0)) / 2


@dataclass(frozen=True)
class VerticalConstruction:
    """Materialized vertical construction (mirrors the imperative trait arrays
    and the ``construction/golden.json`` contract). Optional metric fields are
    ``None`` when the coordinate kind does not define them."""

    levels: np.ndarray
    centers: np.ndarray
    widths: np.ndarray
    cell_volume: np.ndarray
    metric_dz: np.ndarray
    metric_z: np.ndarray
    metric_sigma: np.ndarray | None
    metric_pressure: np.ndarray | None
    metric_ak: np.ndarray | None
    metric_bk: np.ndarray | None
    neighbor_minus: np.ndarray
    neighbor_plus: np.ndarray
    boundary_lower: np.ndarray
    boundary_upper: np.ndarray
    ak: np.ndarray
    bk: np.ndarray
    p0: float


def _faq_levels(grid, np_dtype: type) -> np.ndarray:
    """Synthesize the interface coordinates via the M1 FAQ, branching on the
    coordinate kind exactly as the imperative generator does."""
    coord = grid.coordinate
    nz = int(grid.nz)
    if coord == "eta":
        p0 = float(grid.p0)
        out = np.empty(nz + 1, dtype=np_dtype)
        for k in range(nz + 1):
            out[k] = eval_coeff(
                _ETA, {"ak": float(grid.ak[k]), "p0": p0, "bk": float(grid.bk[k])}
            )
        return out
    if coord in ("sigma", "hybrid_sigma_theta"):
        cand = np.empty(nz + 1, dtype=np_dtype)
        for k in range(nz + 1):
            cand[k] = eval_coeff(_SIGMA, {"k": float(k), "nz": float(nz)})
        levels = np.asarray(grid.levels, dtype=np_dtype)
        return cand if np.array_equal(cand, levels) else levels.copy()
    # z / theta / z_star — supplied verbatim as DATA
    return np.asarray(grid.levels, dtype=np_dtype).copy()


def _faq_centers_widths(levels: np.ndarray, np_dtype: type) -> tuple[np.ndarray, np.ndarray]:
    """Elementwise midpoint / absolute-difference over the ``layers`` interval."""
    n = len(levels) - 1
    centers = np.empty(n, dtype=np_dtype)
    widths = np.empty(n, dtype=np_dtype)
    for k in range(n):
        b = {"a": float(levels[k]), "b": float(levels[k + 1])}
        centers[k] = eval_coeff(_AVG2, b)
        widths[k] = eval_coeff(_WIDTH, b)
    return centers, widths


def _faq_layer_avg(coeff, n: int, np_dtype: type) -> np.ndarray:
    """Per-layer two-point average of an interface-defined coefficient."""
    out = np.empty(n, dtype=np_dtype)
    for k in range(n):
        out[k] = eval_coeff(_AVG2, {"a": float(coeff[k]), "b": float(coeff[k + 1])})
    return out


def vertical_construction_faq(grid) -> VerticalConstruction:
    """Materialize a ``VerticalGrid``'s construction from the declarative FAQ,
    routing all arithmetic through the landed ESS AST evaluator (``eval_coeff``).
    """
    np_dtype = np.dtype(grid.dtype).type
    levels = _faq_levels(grid, np_dtype)
    centers, widths = _faq_centers_widths(levels, np_dtype)
    n = len(widths)

    cell_volume = widths.copy()
    metric_dz = widths.copy()
    metric_z = centers.copy()

    sigma_like = grid.coordinate in ("sigma", "hybrid_sigma_theta", "eta")
    metric_sigma = centers.copy() if sigma_like else None

    has_ak = np.asarray(grid.ak).size > 0
    has_bk = np.asarray(grid.bk).size > 0
    p0 = float(grid.p0)

    metric_pressure: np.ndarray | None
    if has_ak and has_bk:
        metric_pressure = np.empty(n, dtype=np_dtype)
        for k in range(n):
            metric_pressure[k] = eval_coeff(
                _PRESSURE,
                {
                    "ak_lo": float(grid.ak[k]),
                    "bk_lo": float(grid.bk[k]),
                    "ak_hi": float(grid.ak[k + 1]),
                    "bk_hi": float(grid.bk[k + 1]),
                    "p0": p0,
                },
            )
    else:
        metric_pressure = None
    metric_ak = _faq_layer_avg(grid.ak, n, np_dtype) if has_ak else None
    metric_bk = _faq_layer_avg(grid.bk, n, np_dtype) if has_bk else None

    neighbor_minus = np.array([k - 1 if k > 0 else -1 for k in range(n)], dtype=np.int64)
    neighbor_plus = np.array([k + 1 if k + 1 < n else -1 for k in range(n)], dtype=np.int64)
    boundary_lower = np.array([k == 0 for k in range(n)], dtype=bool)
    boundary_upper = np.array([k == n - 1 for k in range(n)], dtype=bool)

    return VerticalConstruction(
        levels=levels,
        centers=centers,
        widths=widths,
        cell_volume=cell_volume,
        metric_dz=metric_dz,
        metric_z=metric_z,
        metric_sigma=metric_sigma,
        metric_pressure=metric_pressure,
        metric_ak=metric_ak,
        metric_bk=metric_bk,
        neighbor_minus=neighbor_minus,
        neighbor_plus=neighbor_plus,
        boundary_lower=boundary_lower,
        boundary_upper=boundary_upper,
        ak=np.asarray(grid.ak, dtype=np_dtype).copy(),
        bk=np.asarray(grid.bk, dtype=np_dtype).copy(),
        p0=float(grid.p0),
    )
