"""DUO icosahedral triangular grid accessor runtime (Heikes et al. 2023).

Conforms to ``docs/GRIDS_API.md`` §2.4 (Python signature), §3.2 (return type),
§7 (common fields), and §10 (loader contract).

Subdivision level ``r`` yields ``20 * 4**r`` triangular cells,
``10 * 4**r + 2`` vertices, and ``30 * 4**r`` edges. All vertices lie on the
sphere of radius ``R``.

Per the 2026-04-20 scope correction the ``.esm`` lowering is a small declarative
config (family + level + loader ref). The geometry is materialized at runtime
inside the generator call (loader-backed families are eagerly resolved per
§10) and then exposed via accessors:

* :meth:`DuoGrid.cell_centers` - cell-center ``(lon, lat)``
* :meth:`DuoGrid.neighbors` - triangle-edge neighbor indices (0-based)
* :meth:`DuoGrid.metric_eval` - per-cell scalar metrics (area, lon/lat, x/y/z)

Only ``builtin://icosahedral/<level>`` loader paths are recognized today; the
``.duo`` file reader lands with the EarthSciSerialization file-format spec
(future bead).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .duo_primal_geometry_faq import primal_geometry_faq
from .duo_subdivide_faq import duo_subdivide_faq
from .duo_topology_faq import primal_topology_faq

__all__ = ["duo", "DuoGrid", "DuoLoader"]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}
_DUO_FAMILY_VERSION = "1.0.0"


@dataclass(frozen=True)
class DuoLoader:
    """Loader spec per ``docs/GRIDS_API.md`` §10.

    ``path`` may be a filesystem path or a ``builtin://icosahedral/<level>``
    URI. ``reader`` selects the decoder; ``check`` controls strictness.
    """

    path: str
    reader: str = "auto"
    check: str = "strict"


def _coerce_loader(loader) -> DuoLoader:
    if isinstance(loader, DuoLoader):
        return loader
    if isinstance(loader, Mapping):
        if "path" not in loader:
            raise TypeError("duo: loader.path is required")
        return DuoLoader(
            path=str(loader["path"]),
            reader=str(loader.get("reader", "auto")),
            check=str(loader.get("check", "strict")),
        )
    raise TypeError(
        f"duo: loader must be DuoLoader or mapping with 'path', got {type(loader).__name__}"
    )


def _parse_builtin_level(path: str) -> int | None:
    prefix = "builtin://icosahedral/"
    if not path.startswith(prefix):
        return None
    tail = path[len(prefix):]
    try:
        lvl = int(tail)
    except ValueError as exc:
        raise ValueError(
            f"duo: cannot parse level from loader path {path!r}"
        ) from exc
    if lvl < 0:
        raise ValueError(f"duo: subdivision level must be >= 0, got {lvl}")
    return lvl


def _resolve_loader_level(loader: DuoLoader) -> int:
    lvl = _parse_builtin_level(loader.path)
    if lvl is not None:
        return lvl
    if loader.reader in ("duo_mesh", "auto"):
        raise ValueError(
            "duo: .duo mesh-file reader not yet implemented - pending "
            "EarthSciSerialization file-format spec. Use "
            "builtin://icosahedral/<level> in the meantime."
        )
    raise ValueError(
        f"duo: unrecognized loader path {loader.path!r} with reader={loader.reader!r}"
    )


# --- Base icosahedron --------------------------------------------------------


def _icosahedron_faces() -> np.ndarray:
    """Return the 20 triangular faces as ``(3, 20)`` of 0-based vertex indices.

    Winding is outward (right-hand rule gives outward normal, so spherical
    triangle area is positive). Ordering mirrors the Julia binding shifted
    from 1-based to 0-based indices.
    """
    faces_1based = np.array(
        [
            [1, 9, 2],
            [1, 2, 11],
            [1, 11, 6],
            [1, 6, 5],
            [1, 5, 9],
            [9, 5, 10],
            [9, 10, 7],
            [9, 7, 2],
            [2, 7, 8],
            [2, 8, 11],
            [11, 8, 12],
            [11, 12, 6],
            [6, 12, 3],
            [6, 3, 5],
            [5, 3, 10],
            [10, 3, 4],
            [10, 4, 7],
            [7, 4, 8],
            [8, 4, 12],
            [12, 4, 3],
        ],
        dtype=np.int64,
    )
    return (faces_1based - 1).T.copy()  # shape (3, 20)


# --- Grid class --------------------------------------------------------------


_METRIC_NAMES = frozenset({"area", "lon", "lat", "x", "y", "z"})


class DuoGrid:
    """DUO icosahedral triangular grid (see module docstring).

    Attributes
    ----------
    level : int
        Subdivision level ``r`` (``r=0`` is the bare icosahedron).
    R : float
        Sphere radius in the configured units.
    dtype : str
        ``"float64"`` or ``"float32"`` (per §2.4).
    ghosts : int
        Halo cell width (§2.7). Always ``0`` for unstructured meshes today.
    vertices : np.ndarray
        ``(3, n_vertices)`` cartesian coordinates scaled by ``R``.
    faces : np.ndarray
        ``(3, n_cells)`` 0-based vertex indices per triangular cell.
    lon, lat : np.ndarray
        ``(n_cells,)`` cell-center geographic coords (radians).
    cell_cart : np.ndarray
        ``(3, n_cells)`` cell-center cartesian coordinates on the sphere.
    area : np.ndarray
        ``(n_cells,)`` spherical-triangle areas (m^2 if ``R`` is in meters).
    edges : np.ndarray
        ``(2, n_edges)`` sorted vertex-pair per edge.
    cell_neighbors : np.ndarray
        ``(3, n_cells)`` neighbor cell indices across each edge (``-1`` if
        boundary; a closed icosahedral mesh has no boundary).
    vertex_faces : list[list[int]]
        Ragged: face indices incident on each vertex (sorted).
    loader : DuoLoader
        Loader the grid was generated from.
    """

    family = "duo"
    topology = "unstructured"

    def __init__(
        self,
        *,
        level: int,
        R: float,
        dtype: str,
        ghosts: int,
        vertices: np.ndarray,
        faces: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        cell_cart: np.ndarray,
        area: np.ndarray,
        edges: np.ndarray,
        cell_neighbors: np.ndarray,
        vertex_faces: list[list[int]],
        loader: DuoLoader,
    ) -> None:
        self.level = level
        self.R = R
        self.dtype = dtype
        self.ghosts = ghosts
        self.vertices = vertices
        self.faces = faces
        self.lon = lon
        self.lat = lat
        self.cell_cart = cell_cart
        self.area = area
        self.edges = edges
        self.cell_neighbors = cell_neighbors
        self.vertex_faces = vertex_faces
        self.loader = loader

    @property
    def n_cells(self) -> int:
        return int(self.faces.shape[1])

    @property
    def n_vertices(self) -> int:
        return int(self.vertices.shape[1])

    @property
    def n_edges(self) -> int:
        return int(self.edges.shape[1])

    def total_area(self) -> float:
        return float(self.area.sum())

    @property
    def provenance(self) -> dict:
        import earthsci_discretizations

        return {
            "binding": "python",
            "binding_version": earthsci_discretizations.__version__,
            "source": "earthsci_discretizations.grids.duo",
            "family": "duo",
            "version": _DUO_FAMILY_VERSION,
            "level": int(self.level),
            "reader": self.loader.reader,
            "path": self.loader.path,
            "check": self.loader.check,
            "dtype": self.dtype,
        }

    def cell_centers(
        self, c: int | None = None
    ) -> tuple[np.ndarray, np.ndarray] | tuple[float, float]:
        """Cell-center geographic coords.

        With no argument, returns ``(lon, lat)`` arrays of shape ``(n_cells,)``.
        With ``c``, returns the scalar ``(lon, lat)`` pair for that cell.
        """
        if c is None:
            return self.lon, self.lat
        self._check_cell(c)
        return float(self.lon[c]), float(self.lat[c])

    def neighbors(self, c: int) -> tuple[int, int, int]:
        """Return the three cell indices sharing an edge with cell ``c``.

        ``-1`` means the edge has no neighbor; a closed icosahedral mesh has
        none.
        """
        self._check_cell(c)
        return (
            int(self.cell_neighbors[0, c]),
            int(self.cell_neighbors[1, c]),
            int(self.cell_neighbors[2, c]),
        )

    def metric_eval(self, name: str, c: int) -> float:
        """Evaluate a per-cell scalar metric.

        Valid ``name`` values:

        * ``"area"`` - spherical-triangle area
        * ``"lon"``, ``"lat"`` - cell-center geographic coords
        * ``"x"``, ``"y"``, ``"z"`` - cell-center cartesian coords
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"duo: metric_eval: unknown metric {name!r}")
        self._check_cell(c)
        if name == "area":
            return float(self.area[c])
        if name == "lon":
            return float(self.lon[c])
        if name == "lat":
            return float(self.lat[c])
        if name == "x":
            return float(self.cell_cart[0, c])
        if name == "y":
            return float(self.cell_cart[1, c])
        return float(self.cell_cart[2, c])

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering per the 2026-04-20 scope correction.

        Returns a §6-schema-shaped config (family + level + loader ref +
        counts + provenance). No inline geometry arrays.
        """
        return {
            "family": self.family,
            "topology": self.topology,
            "dtype": self.dtype,
            "ghosts": int(self.ghosts),
            "n_cells": self.n_cells,
            "n_vertices": self.n_vertices,
            "n_edges": self.n_edges,
            "options": {
                "R": float(self.R),
                "level": int(self.level),
                "loader": {
                    "path": self.loader.path,
                    "reader": self.loader.reader,
                    "check": self.loader.check,
                },
            },
            "provenance": self.provenance,
            "schema_version": _DUO_FAMILY_VERSION,
        }

    def __repr__(self) -> str:
        return (
            f"DuoGrid(level={self.level}, R={self.R}, "
            f"dtype={self.dtype!r}, n_cells={self.n_cells})"
        )

    # internals -------------------------------------------------------------

    def _check_cell(self, c) -> None:
        if c is None:
            raise TypeError("duo: cell index c required")
        ic = int(c)
        if not 0 <= ic < self.n_cells:
            raise ValueError(f"duo: c out of range [0, {self.n_cells}): {c!r}")


# --- Generator ---------------------------------------------------------------


def duo(
    *,
    loader,
    R: float = 6.371e6,
    dtype="float64",
    ghosts: int = 0,
) -> DuoGrid:
    """Generate a DUO icosahedral triangular grid.

    See ``docs/GRIDS_API.md`` §2.4 (Python signature), §3.2 (return type),
    §7 (common fields), §10 (loader contract). All options are keyword-only.

    Parameters
    ----------
    loader : DuoLoader or mapping
        Loader spec per §10. Must provide ``path``; ``reader`` and ``check``
        default to ``"auto"`` and ``"strict"``. Only
        ``builtin://icosahedral/<level>`` paths are accepted today.
    R : float, optional
        Sphere radius. Default ``6.371e6`` (Earth).
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision for the geometry arrays. Default ``"float64"``.
    ghosts : int, optional
        Halo cell width. Must be ``>= 0``; default ``0``.
    """
    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            dtype = "float64"
        elif dtype == np.float32:
            dtype = "float32"
        else:
            raise ValueError(f"duo: unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(f"duo: dtype must be 'float64' or 'float32', got {dtype!r}")

    if isinstance(ghosts, bool) or not isinstance(ghosts, (int, np.integer)):
        raise TypeError(f"duo: ghosts must be int, got {type(ghosts).__name__}")
    ghosts = int(ghosts)
    if ghosts < 0:
        raise ValueError(f"duo: ghosts must be >= 0, got {ghosts}")

    R = float(R)
    if not (R > 0 and math.isfinite(R)):
        raise ValueError(f"duo: R must be a positive finite number, got {R!r}")

    ldr = _coerce_loader(loader)
    level = _resolve_loader_level(ldr)

    np_dtype = _DTYPE_MAP[dtype]

    # ---- Declarative construction via the ESS front-door (esd-un6 / W2) --------
    # The mesh is materialized through the landed ``earthsci_toolkit`` evaluators:
    #   * D3 subdivision  -> ``duo_subdivide_faq`` (every coordinate via
    #     ``eval_coeff``; vertex/face numbering kept in the imperative encounter
    #     order required for the primal-vertex-indexed Voronoi-dual byte-identity),
    #   * D2a primal geometry -> ``primal_geometry_faq`` (``eval_coeff``),
    #   * D1a primal topology -> ``primal_topology_faq`` (value-invention through the
    #     relational skolem/distinct/rank/equijoin primitives; edges in the ESS
    #     canonical sorted order).
    # Edges adopt the canonical (sorted) numbering (no test pins edge order; mirrors
    # the Julia W1 path, ``src/grids/duo.jl``).
    V_unit, F = duo_subdivide_faq(level)  # D3: unit-sphere vertices + faces
    Nv = V_unit.shape[1]

    pg = primal_geometry_faq(V_unit, F, R, np_dtype)  # D2a: cell_cart, lon, lat, area
    cell_cart = pg["cell_cart"]
    lon = pg["lon"]
    lat = pg["lat"]
    area = pg["area"]

    # Scale vertex array to radius R for downstream consumers.
    vertices = np.empty((3, Nv), dtype=np_dtype)
    vertices[0, :] = R * V_unit[0, :]
    vertices[1, :] = R * V_unit[1, :]
    vertices[2, :] = R * V_unit[2, :]

    # D1a value-invention topology through the ESS front-door: edge set + dense
    # numbering (canonical sorted order), cell_neighbors, and per-vertex incident
    # faces. cell_neighbors (cell×slot) and vertex_faces are byte-identical
    # regardless of edge numbering.
    topo = primal_topology_faq(F, Nv)
    edges = topo["edges"]
    cell_neighbors = topo["cell_neighbors"]
    vf = topo["vertex_faces"]

    return DuoGrid(
        level=level,
        R=R,
        dtype=dtype,
        ghosts=ghosts,
        vertices=vertices,
        faces=F,
        lon=lon,
        lat=lat,
        cell_cart=cell_cart,
        area=area,
        edges=edges,
        cell_neighbors=cell_neighbors,
        vertex_faces=vf,
        loader=ldr,
    )
