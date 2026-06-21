"""Arakawa staggered-grid construction via the landed M1 elementwise FAQ path
(``eval_coeff``, :mod:`earthsci_discretizations.rules`).

Python mirror of ``src/arakawa_faq.jl`` (esd-3we, S4). Declarative companion:
``discretizations/grids/arakawa/rules/arakawa_construction.esm``. The construction
is a semiring-FAQ (RFC ``semiring-faq-unified-ir`` §5.1/§5.2: cartesian interval
index-sets + elementwise ``arrayop`` bodies — the M1 path, NO value-invention).

Arakawa staggering over a cartesian base is the cartesian product of two staggered
1-D affine maps per axis: a **center** map ``center(i) = lo + (i-0.5)*d`` and a
**node** map ``node(i) = lo + (i-1)*d``. The four staggered locations are pure
structural gathers (row-major) of those two maps. This module routes every piece
of ARITHMETIC — ``dx``/``dy``, the two axis maps, the cell-volume product —
through ESS's single AST pathway (``eval_coeff``); the structural arrays (the
cartesian-product gather, the row-major neighbor linearization with the ``-1``
off-grid sentinel, the boundary masks, and the static A–E variable-location
table) are pure index logic. Byte-identical to the imperative accessors and to the
committed Julia-reference
``tests/conformance/grids/arakawa/construction/golden.json``.

Takes ``(base, stagger)`` (a ``CartesianBase`` and a stagger label) rather than an
``ArakawaGrid``; the geometry/neighbors/boundary are stagger-independent and the
variable-location table is the only stagger-dependent output.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..rules import eval_coeff

__all__ = ["ArakawaConstruction", "arakawa_construction_faq"]


def _mk(op: str, *args) -> dict:
    return {"op": op, "args": list(args)}


_D = _mk("/", _mk("-", "hi", "lo"), "n")  # (hi - lo) / n
_NODE = _mk("+", "lo", _mk("*", _mk("-", "i", 1), "d"))  # lo + (i - 1)*d
_CENTER = _mk("+", "lo", _mk("*", _mk("-", "i", 0.5), "d"))  # lo + (i - 0.5)*d
_PROD = _mk("*", "a", "b")  # a * b


@dataclass(frozen=True)
class ArakawaConstruction:
    """Materialized arakawa construction (mirrors the imperative accessors and the
    ``construction/golden.json`` contract). Per-cell arrays are over the
    ``(nx, ny)`` cell-center set (row-major); each location's coordinate pair is
    flattened row-major over that location's shape. Location names are the golden
    strings ``cell_center`` / ``u_edge`` / ``v_edge`` / ``corner``."""

    nx: int
    ny: int
    dx: float
    dy: float
    cell_volume: float
    coords_cell_center: tuple[np.ndarray, np.ndarray]
    coords_u_edge: tuple[np.ndarray, np.ndarray]
    coords_v_edge: tuple[np.ndarray, np.ndarray]
    coords_corner: tuple[np.ndarray, np.ndarray]
    neighbor_x_minus: np.ndarray
    neighbor_x_plus: np.ndarray
    neighbor_y_minus: np.ndarray
    neighbor_y_plus: np.ndarray
    boundary_x_lower: np.ndarray
    boundary_x_upper: np.ndarray
    boundary_y_lower: np.ndarray
    boundary_y_upper: np.ndarray
    stagger: str
    rotated: bool
    variable_locations: tuple[str, str, str]


def _center_map(n: int, lo: float, d: float, np_dtype: type) -> np.ndarray:
    out = np.empty(n, dtype=np_dtype)
    for i in range(n):
        out[i] = eval_coeff(_CENTER, {"lo": lo, "i": float(i + 1), "d": d})
    return out


def _node_map(np1: int, lo: float, d: float, np_dtype: type) -> np.ndarray:
    out = np.empty(np1, dtype=np_dtype)
    for i in range(np1):
        out[i] = eval_coeff(_NODE, {"lo": lo, "i": float(i + 1), "d": d})
    return out


def _gather(
    xmap: np.ndarray, ymap: np.ndarray, ni: int, nj: int, np_dtype: type
) -> tuple[np.ndarray, np.ndarray]:
    """Row-major gather: ``out[j*ni + i] = (xmap[i], ymap[j])``."""
    xs = np.tile(np.asarray(xmap[:ni]), nj).astype(np_dtype)
    ys = np.repeat(np.asarray(ymap[:nj]), ni).astype(np_dtype)
    return xs, ys


def _axis_neighbor(nx: int, ny: int, axis_x: bool, offset: int) -> np.ndarray:
    out = np.empty(nx * ny, dtype=np.int64)
    for j in range(ny):
        for i in range(nx):
            ii = i + offset if axis_x else i
            jj = j if axis_x else j + offset
            out[j * nx + i] = (jj * nx + ii) if (0 <= ii < nx and 0 <= jj < ny) else -1
    return out


def _axis_boundary(nx: int, ny: int, axis_x: bool, lower: bool) -> np.ndarray:
    target = 0 if lower else (nx - 1 if axis_x else ny - 1)
    out = np.zeros(nx * ny, dtype=bool)
    for j in range(ny):
        for i in range(nx):
            v = i if axis_x else j
            out[j * nx + i] = v == target
    return out


# Static A–E variable-location table, reconstructed from the stagger label
# INDEPENDENT of the imperative table so the conformance test is a genuine
# cross-check. h is always cell_center.
_VARIABLE_LOCATIONS_UV = {
    "A": ("cell_center", "cell_center"),
    "B": ("corner", "corner"),
    "C": ("u_edge", "v_edge"),
    "D": ("v_edge", "u_edge"),
    "E": ("corner", "corner"),
}


def arakawa_construction_faq(base, stagger: str) -> ArakawaConstruction:
    """Materialize the arakawa construction (over a ``CartesianBase``) from the
    declarative FAQ, routing all arithmetic through ``eval_coeff``."""
    np_dtype = np.float64
    nx = int(base.nx)
    ny = int(base.ny)
    xlo, xhi = float(base.xlo), float(base.xhi)
    ylo, yhi = float(base.ylo), float(base.yhi)

    dx = eval_coeff(_D, {"hi": xhi, "lo": xlo, "n": float(nx)})
    dy = eval_coeff(_D, {"hi": yhi, "lo": ylo, "n": float(ny)})

    center_x = _center_map(nx, xlo, dx, np_dtype)
    center_y = _center_map(ny, ylo, dy, np_dtype)
    node_x = _node_map(nx + 1, xlo, dx, np_dtype)
    node_y = _node_map(ny + 1, ylo, dy, np_dtype)

    coords_cell_center = _gather(center_x, center_y, nx, ny, np_dtype)
    coords_u_edge = _gather(node_x, center_y, nx + 1, ny, np_dtype)
    coords_v_edge = _gather(center_x, node_y, nx, ny + 1, np_dtype)
    coords_corner = _gather(node_x, node_y, nx + 1, ny + 1, np_dtype)

    cell_volume = eval_coeff(_PROD, {"a": dx, "b": dy})

    u_loc, v_loc = _VARIABLE_LOCATIONS_UV[stagger]

    return ArakawaConstruction(
        nx=nx,
        ny=ny,
        dx=dx,
        dy=dy,
        cell_volume=cell_volume,
        coords_cell_center=coords_cell_center,
        coords_u_edge=coords_u_edge,
        coords_v_edge=coords_v_edge,
        coords_corner=coords_corner,
        neighbor_x_minus=_axis_neighbor(nx, ny, True, -1),
        neighbor_x_plus=_axis_neighbor(nx, ny, True, +1),
        neighbor_y_minus=_axis_neighbor(nx, ny, False, -1),
        neighbor_y_plus=_axis_neighbor(nx, ny, False, +1),
        boundary_x_lower=_axis_boundary(nx, ny, True, True),
        boundary_x_upper=_axis_boundary(nx, ny, True, False),
        boundary_y_lower=_axis_boundary(nx, ny, False, True),
        boundary_y_upper=_axis_boundary(nx, ny, False, False),
        stagger=stagger,
        rotated=(stagger == "E"),
        variable_locations=("cell_center", u_loc, v_loc),
    )
