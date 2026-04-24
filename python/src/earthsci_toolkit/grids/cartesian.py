"""Cartesian grid accessor runtime (1D/2D/3D, uniform + non-uniform).

Conforms to the cross-binding contract in ``docs/GRIDS_API.md`` §2.4, §3.2, §7.

Per the 2026-04-20 scope correction, the ``.esm`` lowering is a small
declarative config (family + dimensions + extents, plus explicit edges for
non-uniform axes), NOT a serialized geometry blob. Geometry is derived on
demand by the accessors:

* :meth:`CartesianGrid.cell_centers` - cell-center coordinates
* :meth:`CartesianGrid.cell_widths` - per-axis cell widths
* :meth:`CartesianGrid.cell_volume` - cell measure (length / area / volume)
* :meth:`CartesianGrid.neighbors` - axis-aligned neighbor indices
* :meth:`CartesianGrid.metric_eval` - named metric fields (volume, dx/dy/dz,
  face_area_*, jacobian, g)

The payload shape (field names, order, types) mirrors the Julia reference
binding so that cross-binding `.esm` round-trips compare byte-identical
after canonicalization, modulo the per-binding ``provenance`` block.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = ["cartesian", "CartesianGrid"]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}

_AXIS_NAMES = ("x", "y", "z")

_METRIC_NAMES = frozenset({
    "volume",
    "jacobian",
    "g",
    "dx", "dy", "dz",
    "face_area_x", "face_area_y", "face_area_z",
})


def _is_int(value: Any) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, np.integer))


def _coerce_int(name: str, value: Any) -> int:
    if not _is_int(value):
        raise TypeError(f"cartesian: {name} must be int, got {type(value).__name__}")
    return int(value)


def _coerce_finite_float(name: str, value: Any) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"cartesian: {name} must be a finite number, got {value!r}") from exc
    if not math.isfinite(v):
        raise ValueError(f"cartesian: {name} must be finite, got {value!r}")
    return v


def _coerce_dtype(dtype: Any) -> str:
    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            return "float64"
        if dtype == np.float32:
            return "float32"
        raise ValueError(f"cartesian: unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(
            f"cartesian: dtype must be 'float64' or 'float32', got {dtype!r}"
        )
    return dtype


def _normalize_extent(extent: Any, ndim: int) -> tuple[tuple[float, float], ...]:
    """Normalize ``extent`` into ``((lo, hi), ...)`` of length ``ndim``."""
    if isinstance(extent, (str, bytes)):
        raise TypeError(
            f"cartesian: extent must be a Matrix or Sequence of (lo, hi); "
            f"got {type(extent).__name__}"
        )
    # 2D array form: shape (2, ndim)
    if isinstance(extent, np.ndarray) and extent.ndim == 2:
        if extent.shape[0] != 2:
            raise ValueError(
                f"cartesian: extent matrix must have shape (2, ndim); got {extent.shape}"
            )
        if extent.shape[1] != ndim:
            raise ValueError(
                f"cartesian: extent matrix has {extent.shape[1]} cols, expected {ndim}"
            )
        return tuple(
            (float(extent[0, d]), float(extent[1, d])) for d in range(ndim)
        )
    # Sequence form: length N of (lo, hi) pairs — possibly a nested list/tuple
    # that looks like a 2xN matrix.
    if not isinstance(extent, Sequence) and not isinstance(extent, np.ndarray):
        raise TypeError(
            f"cartesian: extent must be a Matrix or Sequence of (lo, hi); "
            f"got {type(extent).__name__}"
        )
    extent_list = list(extent)
    # Detect 2xN list/tuple form: top-level length 2, each row length ndim.
    if (
        len(extent_list) == 2
        and ndim != 2
        and all(isinstance(row, (list, tuple, np.ndarray)) for row in extent_list)
        and all(len(row) == ndim for row in extent_list)
    ):
        lo_row, hi_row = extent_list
        return tuple(
            (float(lo_row[d]), float(hi_row[d])) for d in range(ndim)
        )
    # Ambiguity for ndim == 2: [[lo1, lo2], [hi1, hi2]] vs [(lo1, hi1), (lo2, hi2)].
    # Disambiguate: in 2D sequence form, each entry is a 2-element (lo, hi) pair,
    # which is true for BOTH interpretations. Prefer the vector-of-pairs reading
    # (which is the primary Julia form per docs/GRIDS_API.md §2.7 extent row).
    if len(extent_list) != ndim:
        raise ValueError(
            f"cartesian: extent sequence has length {len(extent_list)}, expected {ndim}"
        )
    out: list[tuple[float, float]] = []
    for d, pair in enumerate(extent_list):
        if not isinstance(pair, (list, tuple, np.ndarray)) or len(pair) != 2:
            raise ValueError(
                f"cartesian: extent[{d}] must be a 2-element (lo, hi); got {pair!r}"
            )
        out.append((float(pair[0]), float(pair[1])))
    return tuple(out)


def _uniform_axis(
    n: int, lo: float, hi: float, dtype: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if n < 1:
        raise ValueError(f"cartesian: cell count per axis must be >= 1; got {n}")
    if not (hi > lo):
        raise ValueError(f"cartesian: extent must satisfy hi > lo; got lo={lo}, hi={hi}")
    np_dtype = _DTYPE_MAP[dtype]
    dx = (hi - lo) / float(n)
    edges = np.asarray([lo + i * dx for i in range(n + 1)], dtype=np_dtype)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    return edges, centers, widths


def _nonuniform_axis(
    edges_in: Sequence[float], dtype: str
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    if isinstance(edges_in, (str, bytes)):
        raise TypeError("cartesian: edges entry must be a sequence of floats")
    edges_arr = np.asarray(list(edges_in), dtype=np.float64)
    if edges_arr.ndim != 1:
        raise ValueError(
            f"cartesian: edges entry must be 1-D, got shape {edges_arr.shape}"
        )
    if edges_arr.shape[0] < 2:
        raise ValueError(
            f"cartesian: non-uniform edges must have >= 2 entries; got {edges_arr.shape[0]}"
        )
    if not np.all(np.isfinite(edges_arr)):
        raise ValueError("cartesian: edges must be finite")
    if not np.all(np.diff(edges_arr) > 0):
        raise ValueError("cartesian: edges must be strictly increasing")
    np_dtype = _DTYPE_MAP[dtype]
    edges = edges_arr.astype(np_dtype, copy=False)
    n = int(edges.shape[0] - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    return n, edges, centers, widths


def _is_uniform(widths: np.ndarray) -> bool:
    if widths.shape[0] <= 1:
        return True
    w0 = float(widths[0])
    tol = max(np.finfo(widths.dtype).eps, np.finfo(widths.dtype).eps * abs(w0)) * 8.0
    return bool(np.all(np.abs(widths - w0) <= tol))


class CartesianGrid:
    """Uniform or non-uniform Cartesian grid in 1D/2D/3D.

    See module docstring and ``docs/GRIDS_API.md`` §2.4 / §3.2 / §7.

    Attributes
    ----------
    n : tuple of int
        Cells per axis, length ``ndim``.
    extent : tuple of (float, float)
        Per-axis ``(lo, hi)`` bounding box, length ``ndim``.
    edges : tuple of np.ndarray
        Per-axis edge coordinates, each of length ``n[d] + 1``.
    centers : tuple of np.ndarray
        Per-axis cell-center coordinates, each of length ``n[d]``.
    widths : tuple of np.ndarray
        Per-axis cell widths, each of length ``n[d]``.
    uniform : tuple of bool
        Per-axis flag indicating uniform spacing (within dtype ULP tolerance).
    ghosts : int
        Halo cell width.
    dtype : str
        Element precision, ``"float64"`` or ``"float32"``.
    """

    family = "cartesian"
    topology = "rectilinear"

    def __init__(
        self,
        *,
        n: tuple[int, ...],
        extent: tuple[tuple[float, float], ...],
        edges: tuple[np.ndarray, ...],
        centers: tuple[np.ndarray, ...],
        widths: tuple[np.ndarray, ...],
        uniform: tuple[bool, ...],
        ghosts: int,
        dtype: str,
    ) -> None:
        self.n = n
        self.extent = extent
        self.edges = edges
        self.centers = centers
        self.widths = widths
        self.uniform = uniform
        self.ghosts = ghosts
        self.dtype = dtype

    @property
    def ndim(self) -> int:
        return len(self.n)

    @property
    def n_cells(self) -> int:
        out = 1
        for ni in self.n:
            out *= int(ni)
        return out

    @property
    def provenance(self) -> dict:
        import earthsci_toolkit

        return {
            "binding": "python",
            "binding_version": earthsci_toolkit.__version__,
            "source": "earthsci_toolkit.grids.cartesian",
        }

    def cell_centers(self, *idx: int) -> tuple[float, ...]:
        """Return the cell-center coordinates at 0-based index ``(i, [j, [k]])``."""
        self._check_idx(idx)
        return tuple(float(self.centers[d][idx[d]]) for d in range(self.ndim))

    def cell_widths(self, *idx: int) -> tuple[float, ...]:
        """Return per-axis cell widths at 0-based index ``(i, [j, [k]])``."""
        self._check_idx(idx)
        return tuple(float(self.widths[d][idx[d]]) for d in range(self.ndim))

    def cell_volume(self, *idx: int) -> float:
        """Return the cell measure: length (1D), area (2D), or volume (3D)."""
        self._check_idx(idx)
        out = 1.0
        for d in range(self.ndim):
            out *= float(self.widths[d][idx[d]])
        return out

    def neighbors(self, *idx: int) -> dict[tuple[int, int], tuple[int, ...]]:
        """Axis-aligned neighbors as ``{(axis, side): neighbor_idx}``.

        ``axis`` is 0-based in ``range(ndim)``; ``side`` is ``-1`` or ``+1``.
        Boundary cells with no neighbor on a given side are omitted (keeping
        the return type concrete, per the Julia reference semantics).
        """
        self._check_idx(idx)
        out: dict[tuple[int, int], tuple[int, ...]] = {}
        for d in range(self.ndim):
            i_d = int(idx[d])
            if i_d > 0:
                out[(d, -1)] = tuple(
                    (i_d - 1) if k == d else int(idx[k])
                    for k in range(self.ndim)
                )
            if i_d + 1 < int(self.n[d]):
                out[(d, +1)] = tuple(
                    (i_d + 1) if k == d else int(idx[k])
                    for k in range(self.ndim)
                )
        return out

    def metric_eval(
        self, name: str, *idx: int
    ) -> float | tuple[tuple[float, ...], ...]:
        """Evaluate a named metric at cell ``(i, [j, [k]])``.

        Valid names:

        * ``"volume"`` - cell measure.
        * ``"jacobian"`` - ``det(g)^(1/2)`` = 1 for Cartesian.
        * ``"g"`` - metric tensor (identity) as a nested tuple of shape
          ``(ndim, ndim)``.
        * ``"dx"``, ``"dy"``, ``"dz"`` - per-axis cell width (axis must exist).
        * ``"face_area_x"``, ``"face_area_y"``, ``"face_area_z"`` - face area
          normal to the named axis. 1D face area = 1; 2D = orthogonal width;
          3D = product of the two orthogonal widths.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"cartesian: unknown metric name: {name!r}")
        self._check_idx(idx)
        if name == "volume":
            return self.cell_volume(*idx)
        if name == "jacobian":
            return 1.0
        if name == "g":
            return tuple(
                tuple(1.0 if i == j else 0.0 for j in range(self.ndim))
                for i in range(self.ndim)
            )
        if name in ("dx", "dy", "dz"):
            d = ("dx", "dy", "dz").index(name)
            self._check_axis(d)
            return float(self.widths[d][idx[d]])
        # face_area_*
        d = ("face_area_x", "face_area_y", "face_area_z").index(name)
        self._check_axis(d)
        return self._face_area(d, idx)

    def to_esm(self) -> dict:
        """Lower the grid to the §6 declarative ``.esm`` form.

        The payload's field set matches the Julia reference binding so that
        ``canonicalize(to_esm_PY(g)) == canonicalize(to_esm_JL(g))``
        holds byte-for-byte at equivalent opts, modulo the ``provenance``
        block (which is binding-specific per GRIDS_API.md §6.4).
        """
        out: dict[str, Any] = {
            "family": self.family,
            "version": "1.0.0",
            "dtype": self.dtype,
            "topology": self.topology,
            "ndim": int(self.ndim),
            "ghosts": int(self.ghosts),
            "n_cells": int(self.n_cells),
            "n": [int(x) for x in self.n],
            "extent": [[float(lo), float(hi)] for (lo, hi) in self.extent],
            "uniform": [bool(u) for u in self.uniform],
            "provenance": self.provenance,
        }
        if not all(self.uniform):
            out["edges"] = [
                [] if self.uniform[d] else [float(e) for e in self.edges[d]]
                for d in range(self.ndim)
            ]
        return out

    def __repr__(self) -> str:
        return (
            f"CartesianGrid(ndim={self.ndim}, n={self.n}, extent={self.extent}, "
            f"dtype={self.dtype!r}, ghosts={self.ghosts})"
        )

    # internals -----------------------------------------------------------

    def _check_axis(self, d: int) -> None:
        if d >= self.ndim:
            raise ValueError(
                f"cartesian: axis {_AXIS_NAMES[d]} not present in {self.ndim}D grid"
            )

    def _check_idx(self, idx: tuple[int, ...]) -> None:
        if len(idx) != self.ndim:
            raise TypeError(
                f"cartesian: expected {self.ndim} indices, got {len(idx)}"
            )
        for d, v in enumerate(idx):
            if not _is_int(v):
                raise TypeError(
                    f"cartesian: index [{d}] must be int, got {type(v).__name__}"
                )
            if not 0 <= int(v) < int(self.n[d]):
                raise ValueError(
                    f"cartesian: index [{d}]={int(v)} out of range [0, {int(self.n[d])})"
                )

    def _face_area(self, axis: int, idx: tuple[int, ...]) -> float:
        a = 1.0
        for d in range(self.ndim):
            if d == axis:
                continue
            a *= float(self.widths[d][idx[d]])
        return a


def cartesian(
    *,
    nx: int | None = None,
    ny: int | None = None,
    nz: int | None = None,
    extent: Any = None,
    edges: Sequence[Sequence[float]] | None = None,
    dtype: Any = "float64",
    ghosts: int = 0,
) -> CartesianGrid:
    """Generate a Cartesian grid (1D/2D/3D, uniform or non-uniform).

    See ``docs/GRIDS_API.md`` §2.4 for the cross-binding signature contract.
    All options are keyword-only.

    Parameters
    ----------
    nx, ny, nz : int, optional
        Cells per axis. ``nx`` is required for the uniform path; ``ny`` / ``nz``
        are required for 2D / 3D respectively. Forbidden when ``edges`` is
        supplied.
    extent : 2xN array or sequence of (lo, hi), optional
        Axis-aligned bounding box. Required for the uniform path. Accepts
        either a ``(2, ndim)`` numpy array / nested list
        (``[[lo...], [hi...]]``) or a length-``ndim`` sequence of ``(lo, hi)``
        pairs. Forbidden when ``edges`` is supplied.
    edges : sequence of 1-D float sequences, optional
        Per-axis explicit edge coordinates for the non-uniform path. Each
        entry is strictly increasing; axis ``d`` has length ``n[d] + 1``.
        Supersedes ``nx/ny/nz`` and ``extent``.
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision. Default ``"float64"``.
    ghosts : int, optional
        Halo cell width. Default ``0``.
    """
    dtype_s = _coerce_dtype(dtype)
    ghosts_i = _coerce_int("ghosts", ghosts)
    if ghosts_i < 0:
        raise ValueError(f"cartesian: ghosts must be >= 0; got {ghosts_i}")

    if edges is not None:
        if extent is not None or nx is not None or ny is not None or nz is not None:
            raise ValueError(
                "cartesian: pass either `edges` (non-uniform) OR `nx`/`ny`/`nz`+"
                "`extent` (uniform), not both"
            )
        if not isinstance(edges, Sequence) or isinstance(edges, (str, bytes)):
            raise TypeError("cartesian: `edges` must be a sequence of edge arrays")
        if len(edges) == 0:
            raise ValueError("cartesian: `edges` must be a non-empty sequence of edge arrays")
        ndim = len(edges)
        if ndim > 3:
            raise ValueError(f"cartesian: ndim must be <= 3; got {ndim}")
        ax_n: list[int] = []
        ax_e: list[np.ndarray] = []
        ax_c: list[np.ndarray] = []
        ax_w: list[np.ndarray] = []
        ax_u: list[bool] = []
        ax_x: list[tuple[float, float]] = []
        for d in range(ndim):
            n_d, e_d, c_d, w_d = _nonuniform_axis(edges[d], dtype_s)
            ax_n.append(n_d)
            ax_e.append(e_d)
            ax_c.append(c_d)
            ax_w.append(w_d)
            ax_u.append(_is_uniform(w_d))
            ax_x.append((float(e_d[0]), float(e_d[-1])))
        return CartesianGrid(
            n=tuple(ax_n),
            extent=tuple(ax_x),
            edges=tuple(ax_e),
            centers=tuple(ax_c),
            widths=tuple(ax_w),
            uniform=tuple(ax_u),
            ghosts=ghosts_i,
            dtype=dtype_s,
        )

    # Uniform path
    if nx is None:
        raise TypeError("cartesian: required option `nx` is missing")
    nx_i = _coerce_int("nx", nx)
    if ny is None and nz is not None:
        raise TypeError("cartesian: required option `ny` is missing for 3D grid")
    if ny is None:
        ndim = 1
        ns: tuple[int, ...] = (nx_i,)
    else:
        ny_i = _coerce_int("ny", ny)
        if nz is None:
            ndim = 2
            ns = (nx_i, ny_i)
        else:
            nz_i = _coerce_int("nz", nz)
            ndim = 3
            ns = (nx_i, ny_i, nz_i)

    if extent is None:
        raise TypeError("cartesian: required option `extent` is missing")
    ext_tuple = _normalize_extent(extent, ndim)

    ax_e2: list[np.ndarray] = []
    ax_c2: list[np.ndarray] = []
    ax_w2: list[np.ndarray] = []
    for d in range(ndim):
        e, c, w = _uniform_axis(ns[d], ext_tuple[d][0], ext_tuple[d][1], dtype_s)
        ax_e2.append(e)
        ax_c2.append(c)
        ax_w2.append(w)

    return CartesianGrid(
        n=ns,
        extent=ext_tuple,
        edges=tuple(ax_e2),
        centers=tuple(ax_c2),
        widths=tuple(ax_w2),
        uniform=tuple(True for _ in range(ndim)),
        ghosts=ghosts_i,
        dtype=dtype_s,
    )
