"""MPAS unstructured Voronoi grid accessor runtime (loader-backed).

Conforms to ``docs/GRIDS_API.md`` §2.4 (Python signature), §3.2 (return type),
§7 (common fields), and §10 (loader contract), and to the mayor's 2026-04-20
scope correction: ``to_esm()`` emits a small declarative config (family,
dimensions, loader ref) — no inline geometry arrays. Geometry is derived on
demand by the accessors (:meth:`MpasGrid.cell_centers`,
:meth:`MpasGrid.neighbors`, :meth:`MpasGrid.cell_area`,
:meth:`MpasGrid.edge_length`, :meth:`MpasGrid.metric_eval`). Cross-binding
conformance is compared at pinned query points.

NetCDF I/O is intentionally not bundled with ``earthsci_discretizations``. Path-based
loading requires the caller to pass a ``reader_fn(path) -> MpasMeshData``;
in-memory construction via :func:`mpas_mesh_data` is the primary path for
tests and host-built meshes. Mirrors the TypeScript binding's ``readerFn``
contract in ``typescript/src/grids/mpas.ts``.

Index convention: cells, edges, and adjacency slots are 0-based. The sentinel
for "no neighbor" / "no edge" in adjacency arrays is ``-1``. This matches the
TypeScript binding.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

__all__ = [
    "mpas",
    "MpasGrid",
    "MpasLoader",
    "MpasMeshData",
    "mpas_mesh_data",
    "check_mesh",
]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}
_MPAS_FAMILY_VERSION = "1.0.0"
_MPAS_SOURCE_SHA = "dsc-uct"
_MPAS_READER_VERSION = "0.1.0"

_VALID_READERS = frozenset({"auto", "nc4", "mpas_mesh"})
_VALID_CHECKS = frozenset({"strict", "lenient"})

_METRIC_CELL = frozenset({"lon", "lat", "area", "x", "y", "z", "n_edges_on_cell"})
_METRIC_EDGE = frozenset({"lon_edge", "lat_edge", "dc_edge", "dv_edge"})
_METRIC_NAMES = _METRIC_CELL | _METRIC_EDGE


@dataclass(frozen=True)
class MpasLoader:
    """Loader spec per ``docs/GRIDS_API.md`` §10.

    Attributes
    ----------
    path : str
        Filesystem path to an MPAS mesh file.
    reader : {"auto", "nc4", "mpas_mesh"}
        Reader selection (default ``"auto"``).
    check : {"strict", "lenient"}
        Mesh validation strictness (default ``"strict"``).
    """

    path: str
    reader: str = "auto"
    check: str = "strict"


def _coerce_loader(loader) -> MpasLoader:
    if isinstance(loader, MpasLoader):
        ldr = loader
    elif isinstance(loader, Mapping):
        if "path" not in loader:
            raise TypeError("mpas: loader.path is required")
        ldr = MpasLoader(
            path=str(loader["path"]),
            reader=str(loader.get("reader", "auto")),
            check=str(loader.get("check", "strict")),
        )
    else:
        raise TypeError(
            f"mpas: loader must be MpasLoader or mapping with 'path', got {type(loader).__name__}"
        )
    if not isinstance(ldr.path, str) or len(ldr.path) == 0:
        raise TypeError("mpas: loader.path must be a non-empty string")
    if ldr.reader not in _VALID_READERS:
        raise ValueError(
            f"mpas: loader.reader must be one of {sorted(_VALID_READERS)}, got {ldr.reader!r}"
        )
    if ldr.check not in _VALID_CHECKS:
        raise ValueError(
            f"mpas: loader.check must be one of {sorted(_VALID_CHECKS)}, got {ldr.check!r}"
        )
    return ldr


# --- Mesh data ---------------------------------------------------------------


@dataclass(frozen=True)
class MpasMeshData:
    """In-memory MPAS Voronoi mesh (validated arrays).

    All arrays use ``float64`` / ``int32`` for cross-binding byte stability
    in the accessor runtime. Adjacency arrays use ``-1`` as the "no neighbor"
    / "no edge" sentinel.
    """

    n_cells: int
    n_edges: int
    n_vertices: int
    max_edges: int
    lon_cell: np.ndarray
    lat_cell: np.ndarray
    x_cell: np.ndarray
    y_cell: np.ndarray
    z_cell: np.ndarray
    area_cell: np.ndarray
    n_edges_on_cell: np.ndarray
    cells_on_cell: np.ndarray  # shape (n_cells, max_edges)
    edges_on_cell: np.ndarray  # shape (n_cells, max_edges)
    lon_edge: np.ndarray
    lat_edge: np.ndarray
    cells_on_edge: np.ndarray  # shape (n_edges, 2)
    dc_edge: np.ndarray
    dv_edge: np.ndarray


def _as_f64_1d(arr, label: str, n: int) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim != 1 or a.shape[0] != n:
        raise ValueError(f"mpas: {label} must be a 1-D length-{n} sequence (got shape {a.shape})")
    return np.ascontiguousarray(a)


def _as_i32_1d(arr, label: str, n: int) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim != 1 or a.shape[0] != n:
        raise ValueError(f"mpas: {label} must be a 1-D length-{n} sequence (got shape {a.shape})")
    if not np.issubdtype(a.dtype, np.integer):
        if not np.all(np.equal(np.mod(a.astype(np.float64), 1.0), 0.0)):
            raise TypeError(f"mpas: {label} entries must be integers")
    return np.ascontiguousarray(a.astype(np.int32))


def _as_i32_2d(arr, label: str, rows: int, cols: int) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim == 1:
        if a.shape[0] != rows * cols:
            raise ValueError(
                f"mpas: {label} must be length {rows * cols} (flat {rows}×{cols}); got {a.shape}"
            )
        a = a.reshape(rows, cols)
    if a.shape != (rows, cols):
        raise ValueError(f"mpas: {label} must be shape ({rows}, {cols}); got {a.shape}")
    if not np.issubdtype(a.dtype, np.integer):
        if not np.all(np.equal(np.mod(a.astype(np.float64), 1.0), 0.0)):
            raise TypeError(f"mpas: {label} entries must be integers")
    return np.ascontiguousarray(a.astype(np.int32))


def mpas_mesh_data(
    *,
    lon_cell,
    lat_cell,
    area_cell,
    n_edges_on_cell,
    cells_on_cell,
    edges_on_cell,
    lon_edge,
    lat_edge,
    cells_on_edge,
    dc_edge,
    dv_edge,
    max_edges: int,
    x_cell=None,
    y_cell=None,
    z_cell=None,
    n_vertices: int = 0,
    R: float = 6.371e6,
) -> MpasMeshData:
    """Validated in-memory MPAS mesh constructor.

    Derives ``x_cell``, ``y_cell``, ``z_cell`` from ``(lon_cell, lat_cell, R)``
    when not supplied.
    """
    if not isinstance(max_edges, (int, np.integer)) or isinstance(max_edges, bool):
        raise TypeError(f"mpas: max_edges must be int, got {type(max_edges).__name__}")
    max_edges = int(max_edges)
    if max_edges <= 0:
        raise ValueError(f"mpas: max_edges must be a positive integer, got {max_edges}")

    lon_cell_arr = np.asarray(lon_cell, dtype=np.float64)
    n_cells = int(lon_cell_arr.shape[0])
    lon_edge_arr = np.asarray(lon_edge, dtype=np.float64)
    n_edges = int(lon_edge_arr.shape[0])

    lon_cell_f = _as_f64_1d(lon_cell, "lon_cell", n_cells)
    lat_cell_f = _as_f64_1d(lat_cell, "lat_cell", n_cells)
    area_cell_f = _as_f64_1d(area_cell, "area_cell", n_cells)
    neoc = _as_i32_1d(n_edges_on_cell, "n_edges_on_cell", n_cells)
    coc = _as_i32_2d(cells_on_cell, "cells_on_cell", n_cells, max_edges)
    eoc = _as_i32_2d(edges_on_cell, "edges_on_cell", n_cells, max_edges)
    lon_edge_f = _as_f64_1d(lon_edge, "lon_edge", n_edges)
    lat_edge_f = _as_f64_1d(lat_edge, "lat_edge", n_edges)
    coe = _as_i32_2d(cells_on_edge, "cells_on_edge", n_edges, 2)
    dc_edge_f = _as_f64_1d(dc_edge, "dc_edge", n_edges)
    dv_edge_f = _as_f64_1d(dv_edge, "dv_edge", n_edges)

    if x_cell is None:
        x_f = R * np.cos(lat_cell_f) * np.cos(lon_cell_f)
    else:
        x_f = _as_f64_1d(x_cell, "x_cell", n_cells)
    if y_cell is None:
        y_f = R * np.cos(lat_cell_f) * np.sin(lon_cell_f)
    else:
        y_f = _as_f64_1d(y_cell, "y_cell", n_cells)
    if z_cell is None:
        z_f = R * np.sin(lat_cell_f)
    else:
        z_f = _as_f64_1d(z_cell, "z_cell", n_cells)

    return MpasMeshData(
        n_cells=n_cells,
        n_edges=n_edges,
        n_vertices=int(n_vertices),
        max_edges=max_edges,
        lon_cell=lon_cell_f,
        lat_cell=lat_cell_f,
        x_cell=x_f,
        y_cell=y_f,
        z_cell=z_f,
        area_cell=area_cell_f,
        n_edges_on_cell=neoc,
        cells_on_cell=coc,
        edges_on_cell=eoc,
        lon_edge=lon_edge_f,
        lat_edge=lat_edge_f,
        cells_on_edge=coe,
        dc_edge=dc_edge_f,
        dv_edge=dv_edge_f,
    )


def check_mesh(mesh: MpasMeshData, strict: bool) -> None:
    """Bounds-check adjacency arrays; in strict mode enforce reciprocity.

    Mirrors ``checkMesh`` in ``typescript/src/grids/mpas.ts``.
    """
    n_cells = mesh.n_cells
    n_edges = mesh.n_edges
    max_edges = mesh.max_edges
    for c in range(n_cells):
        k = int(mesh.n_edges_on_cell[c])
        if k < 0 or k > max_edges:
            raise AssertionError(f"mpas: n_edges_on_cell[{c}]={k} out of [0, {max_edges}]")
        for j in range(k):
            nb = int(mesh.cells_on_cell[c, j])
            if nb < -1 or nb >= n_cells:
                raise AssertionError(
                    f"mpas: cells_on_cell[{c},{j}]={nb} out of [-1, {n_cells - 1}]"
                )
            e = int(mesh.edges_on_cell[c, j])
            if e < -1 or e >= n_edges:
                raise AssertionError(f"mpas: edges_on_cell[{c},{j}]={e} out of [-1, {n_edges - 1}]")
    for e in range(n_edges):
        for s in range(2):
            c = int(mesh.cells_on_edge[e, s])
            if c < -1 or c >= n_cells:
                raise AssertionError(f"mpas: cells_on_edge[{e},{s}]={c} out of [-1, {n_cells - 1}]")
    if not strict:
        return
    for c in range(n_cells):
        k = int(mesh.n_edges_on_cell[c])
        for j in range(k):
            nb = int(mesh.cells_on_cell[c, j])
            if nb < 0:
                continue
            kb = int(mesh.n_edges_on_cell[nb])
            found = False
            for jj in range(kb):
                if int(mesh.cells_on_cell[nb, jj]) == c:
                    found = True
                    break
            if not found:
                raise AssertionError(
                    f"mpas: neighbor symmetry broken: cell {c} -> {nb} but not reverse"
                )


# --- Grid class --------------------------------------------------------------


class MpasGrid:
    """MPAS unstructured Voronoi grid (see module docstring).

    Attributes
    ----------
    family : str
        Always ``"mpas"``.
    topology : str
        Always ``"unstructured"``.
    dtype : str
        ``"float64"`` or ``"float32"`` (per §2.4).
    R : float
        Sphere radius (m).
    ghosts : int
        Halo cell width. Always ``0`` for loader-backed meshes.
    n_cells, n_edges, n_vertices, max_edges : int
        Mesh dimensions.
    mesh : MpasMeshData
        Underlying validated mesh arrays.
    loader : MpasLoader | None
        Loader record if path-based or loader-described; ``None`` for pure
        in-memory construction.
    """

    family = "mpas"
    topology = "unstructured"

    def __init__(
        self,
        *,
        mesh: MpasMeshData,
        R: float,
        dtype: str,
        ghosts: int,
        loader: MpasLoader | None,
    ) -> None:
        self.mesh = mesh
        self.R = R
        self.dtype = dtype
        self.ghosts = ghosts
        self.loader = loader

    @property
    def n_cells(self) -> int:
        return int(self.mesh.n_cells)

    @property
    def n_edges(self) -> int:
        return int(self.mesh.n_edges)

    @property
    def n_vertices(self) -> int:
        return int(self.mesh.n_vertices)

    @property
    def max_edges(self) -> int:
        return int(self.mesh.max_edges)

    @property
    def provenance(self) -> dict:
        import earthsci_discretizations

        ldr = self.loader
        loader_blk = (
            None if ldr is None else {"path": ldr.path, "reader": ldr.reader, "check": ldr.check}
        )
        return {
            "binding": "python",
            "binding_version": earthsci_discretizations.__version__,
            "source": "earthsci_discretizations.grids.mpas",
            "family": "mpas",
            "version": _MPAS_FAMILY_VERSION,
            "source_sha": _MPAS_SOURCE_SHA,
            "reader_version": _MPAS_READER_VERSION,
            "loader": loader_blk,
            "dtype": self.dtype,
        }

    # --- accessors ---------------------------------------------------------

    def cell_centers(self, c: int | None = None):
        """Cell-center geographic coords.

        With no argument, returns ``(lon, lat)`` arrays of shape ``(n_cells,)``.
        With ``c``, returns the scalar ``(lon, lat)`` pair for that cell.
        """
        if c is None:
            return self.mesh.lon_cell, self.mesh.lat_cell
        self._check_cell(c)
        ic = int(c)
        return float(self.mesh.lon_cell[ic]), float(self.mesh.lat_cell[ic])

    def cell_center_cart(self, c: int) -> tuple[float, float, float]:
        """Cell-center cartesian coords (``R``-scaled)."""
        self._check_cell(c)
        ic = int(c)
        return (
            float(self.mesh.x_cell[ic]),
            float(self.mesh.y_cell[ic]),
            float(self.mesh.z_cell[ic]),
        )

    def neighbors(self, c: int) -> list[int]:
        """Indices of cells sharing an edge with cell ``c`` (``-1`` filtered)."""
        self._check_cell(c)
        ic = int(c)
        k = int(self.mesh.n_edges_on_cell[ic])
        out: list[int] = []
        for j in range(k):
            nb = int(self.mesh.cells_on_cell[ic, j])
            if nb >= 0:
                out.append(nb)
        return out

    def cell_area(self, c: int) -> float:
        self._check_cell(c)
        return float(self.mesh.area_cell[int(c)])

    def edge_length(self, e: int) -> float:
        """Voronoi-vertex-to-Voronoi-vertex edge arc length (``dv_edge``)."""
        self._check_edge(e)
        return float(self.mesh.dv_edge[int(e)])

    def cell_distance(self, e: int) -> float:
        """Cell-center-to-cell-center great-circle distance (``dc_edge``)."""
        self._check_edge(e)
        return float(self.mesh.dc_edge[int(e)])

    def total_area(self) -> float:
        return float(self.mesh.area_cell.sum())

    def metric_eval(self, name: str, i: int) -> float:
        """Evaluate a per-cell or per-edge scalar metric.

        Valid cell metrics: ``"lon"``, ``"lat"``, ``"area"``, ``"x"``,
        ``"y"``, ``"z"``, ``"n_edges_on_cell"``. Valid edge metrics:
        ``"lon_edge"``, ``"lat_edge"``, ``"dc_edge"``, ``"dv_edge"``.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"mpas: metric_eval: unknown metric {name!r}")
        if name in _METRIC_CELL:
            self._check_cell(i)
            ii = int(i)
            if name == "lon":
                return float(self.mesh.lon_cell[ii])
            if name == "lat":
                return float(self.mesh.lat_cell[ii])
            if name == "area":
                return float(self.mesh.area_cell[ii])
            if name == "x":
                return float(self.mesh.x_cell[ii])
            if name == "y":
                return float(self.mesh.y_cell[ii])
            if name == "z":
                return float(self.mesh.z_cell[ii])
            return float(self.mesh.n_edges_on_cell[ii])
        self._check_edge(i)
        ii = int(i)
        if name == "lon_edge":
            return float(self.mesh.lon_edge[ii])
        if name == "lat_edge":
            return float(self.mesh.lat_edge[ii])
        if name == "dc_edge":
            return float(self.mesh.dc_edge[ii])
        return float(self.mesh.dv_edge[ii])

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering per the 2026-04-20 scope correction.

        Shape matches the TypeScript binding's ``toESM()``. No inline geometry
        arrays — cross-binding conformance compares accessor outputs at pinned
        query points, not serialized bytes (see GRIDS_API §3.5 and the
        mayor's 2026-04-20 correction on bead dsc-3nw).
        """
        ldr = self.loader
        loader_blk = (
            None if ldr is None else {"path": ldr.path, "reader": ldr.reader, "check": ldr.check}
        )
        return {
            "family": self.family,
            "version": _MPAS_FAMILY_VERSION,
            "dtype": self.dtype,
            "topology": self.topology,
            "ghosts": int(self.ghosts),
            "n_cells": self.n_cells,
            "n_edges": self.n_edges,
            "n_vertices": self.n_vertices,
            "max_edges": self.max_edges,
            "options": {
                "R": float(self.R),
                "loader": loader_blk,
            },
            "provenance": self.provenance,
            "schema_version": _MPAS_FAMILY_VERSION,
        }

    def __repr__(self) -> str:
        return (
            f"MpasGrid(n_cells={self.n_cells}, n_edges={self.n_edges}, "
            f"max_edges={self.max_edges}, R={self.R}, dtype={self.dtype!r})"
        )

    # --- internals ---------------------------------------------------------

    def _check_cell(self, c) -> None:
        if c is None or isinstance(c, bool):
            raise TypeError("mpas: cell index must be int")
        try:
            ic = int(c)
        except (TypeError, ValueError) as exc:
            raise TypeError("mpas: cell index must be int") from exc
        if ic != c and not isinstance(c, np.integer):
            raise TypeError(f"mpas: cell index must be int, got {type(c).__name__}")
        if not 0 <= ic < self.n_cells:
            raise IndexError(f"mpas: invalid cell index {c!r} (expected 0..{self.n_cells - 1})")

    def _check_edge(self, e) -> None:
        if e is None or isinstance(e, bool):
            raise TypeError("mpas: edge index must be int")
        try:
            ie = int(e)
        except (TypeError, ValueError) as exc:
            raise TypeError("mpas: edge index must be int") from exc
        if ie != e and not isinstance(e, np.integer):
            raise TypeError(f"mpas: edge index must be int, got {type(e).__name__}")
        if not 0 <= ie < self.n_edges:
            raise IndexError(f"mpas: invalid edge index {e!r} (expected 0..{self.n_edges - 1})")


# --- Generator ---------------------------------------------------------------


def mpas(
    *,
    mesh: MpasMeshData | None = None,
    loader=None,
    reader_fn=None,
    R: float = 6.371e6,
    dtype="float64",
    ghosts: int = 0,
) -> MpasGrid:
    """Generate an MPAS unstructured Voronoi grid.

    See ``docs/GRIDS_API.md`` §2.4 (Python signature), §3.2 (return type),
    §7 (common fields), §10 (loader contract). All options are keyword-only.

    Parameters
    ----------
    mesh : MpasMeshData, optional
        Pre-loaded in-memory mesh. One of ``mesh`` or ``loader`` is required.
    loader : MpasLoader or mapping, optional
        Loader spec per §10. If ``mesh`` is not supplied, ``reader_fn`` must
        also be provided to translate the loader path into an
        :class:`MpasMeshData`. NetCDF I/O is not bundled — callers provide
        the reader.
    reader_fn : callable, optional
        ``reader_fn(path: str) -> MpasMeshData``. Required when loading from
        a loader path (i.e. ``mesh`` is ``None``).
    R : float, optional
        Sphere radius. Default ``6.371e6`` (Earth).
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision. Default ``"float64"``.
    ghosts : int, optional
        Halo cell width. Must be ``0`` for loader-backed meshes.
    """
    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            dtype = "float64"
        elif dtype == np.float32:
            dtype = "float32"
        else:
            raise ValueError(f"mpas: unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(f"mpas: dtype must be 'float64' or 'float32', got {dtype!r}")

    if isinstance(ghosts, bool) or not isinstance(ghosts, (int, np.integer)):
        raise TypeError(f"mpas: ghosts must be int, got {type(ghosts).__name__}")
    ghosts = int(ghosts)
    if ghosts != 0:
        raise ValueError(f"mpas: ghosts must be 0 for loader-backed grids, got {ghosts}")

    R = float(R)
    if not (R > 0 and math.isfinite(R)):
        raise ValueError(f"mpas: R must be a positive finite number, got {R!r}")

    loader_record: MpasLoader | None = None
    resolved_mesh: MpasMeshData

    if mesh is not None:
        if not isinstance(mesh, MpasMeshData):
            raise TypeError("mpas: 'mesh' must be an MpasMeshData (see mpas_mesh_data())")
        resolved_mesh = mesh
        if loader is not None:
            loader_record = _coerce_loader(loader)
    elif loader is not None:
        loader_record = _coerce_loader(loader)
        if reader_fn is None:
            raise TypeError(
                "mpas: path-based loading requires reader_fn(path) -> MpasMeshData. "
                "NetCDF I/O is not bundled with earthsci_discretizations per "
                "GRIDS_API.md §10; pass reader_fn from your consumer package."
            )
        if not callable(reader_fn):
            raise TypeError("mpas: reader_fn must be callable")
        produced = reader_fn(loader_record.path)
        if not isinstance(produced, MpasMeshData):
            raise TypeError("mpas: reader_fn must return an MpasMeshData instance")
        resolved_mesh = produced
    else:
        raise TypeError(
            "mpas: provide either 'mesh' (MpasMeshData) or 'loader' (MpasLoader) with a reader_fn"
        )

    strict = True if loader_record is None else loader_record.check == "strict"
    check_mesh(resolved_mesh, strict)

    return MpasGrid(
        mesh=resolved_mesh,
        R=R,
        dtype=dtype,
        ghosts=ghosts,
        loader=loader_record,
    )
