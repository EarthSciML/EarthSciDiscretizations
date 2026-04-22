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


def _icosahedron_vertices() -> np.ndarray:
    """Return the 12 canonical icosahedron unit-sphere vertices as ``(3, 12)``.

    Ordering is fixed so subdivision level-0 output is byte-stable across runs
    and platforms and matches the Julia binding (:file:`src/grids/duo.jl`).
    """
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = np.array(
        [
            [0.0, 1.0, phi],
            [0.0, -1.0, phi],
            [0.0, 1.0, -phi],
            [0.0, -1.0, -phi],
            [1.0, phi, 0.0],
            [-1.0, phi, 0.0],
            [1.0, -phi, 0.0],
            [-1.0, -phi, 0.0],
            [phi, 0.0, 1.0],
            [phi, 0.0, -1.0],
            [-phi, 0.0, 1.0],
            [-phi, 0.0, -1.0],
        ],
        dtype=np.float64,
    )
    norms = np.linalg.norm(raw, axis=1, keepdims=True)
    return (raw / norms).T.copy()  # shape (3, 12)


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


# --- Recursive subdivision ---------------------------------------------------


def _subdivide_icosahedron(level: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(vertices (3, Nv), faces (3, Nc))`` after ``level`` subdivisions."""
    V0 = _icosahedron_vertices()
    F0 = _icosahedron_faces()
    verts: list[tuple[float, float, float]] = [
        (float(V0[0, i]), float(V0[1, i]), float(V0[2, i])) for i in range(V0.shape[1])
    ]
    faces: list[tuple[int, int, int]] = [
        (int(F0[0, c]), int(F0[1, c]), int(F0[2, c])) for c in range(F0.shape[1])
    ]

    for _ in range(level):
        cache: dict[tuple[int, int], int] = {}

        def midpoint(a: int, b: int) -> int:
            key = (a, b) if a < b else (b, a)
            idx = cache.get(key)
            if idx is not None:
                return idx
            va = verts[a]
            vb = verts[b]
            mx = va[0] + vb[0]
            my = va[1] + vb[1]
            mz = va[2] + vb[2]
            n = math.sqrt(mx * mx + my * my + mz * mz)
            verts.append((mx / n, my / n, mz / n))
            idx = len(verts) - 1
            cache[key] = idx
            return idx

        new_faces: list[tuple[int, int, int]] = []
        for a, b, c in faces:
            ab = midpoint(a, b)
            bc = midpoint(b, c)
            ca = midpoint(c, a)
            new_faces.append((a, ab, ca))
            new_faces.append((b, bc, ab))
            new_faces.append((c, ca, bc))
            new_faces.append((ab, bc, ca))
        faces = new_faces

    Nv = len(verts)
    Nc = len(faces)
    V = np.empty((3, Nv), dtype=np.float64)
    for i, v in enumerate(verts):
        V[0, i], V[1, i], V[2, i] = v
    F = np.empty((3, Nc), dtype=np.int64)
    for c, f in enumerate(faces):
        F[0, c], F[1, c], F[2, c] = f
    return V, F


# --- Geometry helpers --------------------------------------------------------


def _cart_to_lonlat(x: float, y: float, z: float) -> tuple[float, float]:
    return (math.atan2(y, x), math.asin(max(-1.0, min(1.0, z))))


def _spherical_triangle_area(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    """Spherical excess of the unit-sphere triangle ``(a, b, c)`` via L'Huilier.

    Numerically stable on near-degenerate triangles - matches the Julia
    binding's ``_spherical_triangle_area``.
    """
    da = math.acos(max(-1.0, min(1.0, b[0] * c[0] + b[1] * c[1] + b[2] * c[2])))
    db = math.acos(max(-1.0, min(1.0, c[0] * a[0] + c[1] * a[1] + c[2] * a[2])))
    dc = math.acos(max(-1.0, min(1.0, a[0] * b[0] + a[1] * b[1] + a[2] * b[2])))
    s = 0.5 * (da + db + dc)
    t = (
        math.tan(s / 2)
        * math.tan((s - da) / 2)
        * math.tan((s - db) / 2)
        * math.tan((s - dc) / 2)
    )
    if t < 0.0:
        t = 0.0
    return 4.0 * math.atan(math.sqrt(t))


# --- Connectivity ------------------------------------------------------------


def _build_connectivity(
    faces: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Build ``(edges (2, Ne), cell_neighbors (3, Nc))`` from faces ``(3, Nc)``.

    ``cell_neighbors[k, c]`` is the cell sharing the edge opposite vertex
    ``faces[k, c]`` in cell ``c`` (that edge spans ``faces[(k+1) % 3, c]`` and
    ``faces[(k+2) % 3, c]``). ``-1`` means no neighbor.
    """
    Nc = faces.shape[1]
    edge_map: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for c in range(Nc):
        v1, v2, v3 = int(faces[0, c]), int(faces[1, c]), int(faces[2, c])
        pairs = ((v2, v3), (v3, v1), (v1, v2))
        for k, (a, b) in enumerate(pairs):
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append((c, k))

    Ne = len(edge_map)
    edges = np.empty((2, Ne), dtype=np.int64)
    neighbors = np.full((3, Nc), -1, dtype=np.int64)
    for e, (key, cells) in enumerate(edge_map.items()):
        edges[0, e] = key[0]
        edges[1, e] = key[1]
        if len(cells) == 2:
            (c1, k1), (c2, k2) = cells
            neighbors[k1, c1] = c2
            neighbors[k2, c2] = c1
        elif len(cells) > 2:
            raise AssertionError(
                f"duo: non-manifold edge {key} shared by {len(cells)} cells"
            )
    return edges, neighbors


def _vertex_faces(faces: np.ndarray, Nv: int) -> list[list[int]]:
    vf: list[list[int]] = [[] for _ in range(Nv)]
    Nc = faces.shape[1]
    for c in range(Nc):
        vf[int(faces[0, c])].append(c)
        vf[int(faces[1, c])].append(c)
        vf[int(faces[2, c])].append(c)
    for v in range(Nv):
        vf[v].sort()
    return vf


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
        import earthsci_toolkit

        return {
            "binding": "python",
            "binding_version": earthsci_toolkit.__version__,
            "source": "earthsci_toolkit.grids.duo",
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

    V_unit, F = _subdivide_icosahedron(level)
    Nv = V_unit.shape[1]
    Nc = F.shape[1]

    # Cell centers: normalized centroid of the three unit vertex vectors.
    cell_cart = np.empty((3, Nc), dtype=np_dtype)
    lon = np.empty(Nc, dtype=np_dtype)
    lat = np.empty(Nc, dtype=np_dtype)
    area = np.empty(Nc, dtype=np_dtype)
    R2 = R * R

    for c in range(Nc):
        a_i = int(F[0, c])
        b_i = int(F[1, c])
        c_i = int(F[2, c])
        ax, ay, az = float(V_unit[0, a_i]), float(V_unit[1, a_i]), float(V_unit[2, a_i])
        bx, by, bz = float(V_unit[0, b_i]), float(V_unit[1, b_i]), float(V_unit[2, b_i])
        cx, cy, cz = float(V_unit[0, c_i]), float(V_unit[1, c_i]), float(V_unit[2, c_i])
        mx = ax + bx + cx
        my = ay + by + cy
        mz = az + bz + cz
        n = math.sqrt(mx * mx + my * my + mz * mz)
        ux, uy, uz = mx / n, my / n, mz / n
        cell_cart[0, c] = R * ux
        cell_cart[1, c] = R * uy
        cell_cart[2, c] = R * uz
        lo, la = _cart_to_lonlat(ux, uy, uz)
        lon[c] = lo
        lat[c] = la
        a_unit = _spherical_triangle_area(
            (ax, ay, az), (bx, by, bz), (cx, cy, cz)
        )
        area[c] = a_unit * R2

    vertices = np.empty((3, Nv), dtype=np_dtype)
    vertices[0, :] = R * V_unit[0, :]
    vertices[1, :] = R * V_unit[1, :]
    vertices[2, :] = R * V_unit[2, :]

    edges, cell_neighbors = _build_connectivity(F)
    vf = _vertex_faces(F, Nv)

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
