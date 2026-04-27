#!/usr/bin/env python3
"""
render_doc_plots — generate grid + rule visualizations for the Hugo doc site.

Writes PNG artifacts under ``docs/static/plots/``:

    grids/<family>.png            — typical configuration of each grid family
    rules/<rule>-stencil.png      — stencil + coefficient diagram per rule
    rules/<rule>-convergence.png  — empirical convergence (only for rules
                                    whose Layer-B fixtures are currently
                                    producing — Sec. APPLICABLE below)

Rules whose convergence fixtures depend on in-flight ESS harness extensions
do not produce a convergence plot here; the doc page renders a
``pending`` callout instead.

This is a documentation tool, not a numerics oracle. It uses idealized
mathematical re-implementations of the rule and a manufactured solution to
demonstrate the *shape* of the convergence curve (slope, asymptote). The
authoritative convergence numbers come from the ESS walker; once the in-flight
harness extensions land and the fixtures populate, this script will be
extended to consume the walker output directly.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Layer-B applicability table (mirrors rule-catalog.md "applicable" axis).
# Rules listed here have a working Layer-B fixture today; everything else
# renders a placeholder convergence section.
# ---------------------------------------------------------------------------

APPLICABLE = {
    "centered_2nd_uniform",
    "centered_2nd_uniform_vertical",
    "centered_2nd_uniform_latlon",
    "ppm_reconstruction",
    "upwind_1st",
    "weno5_advection_2d",
}

ALL_RULES = (
    "centered_2nd_uniform",
    "centered_2nd_uniform_vertical",
    "centered_2nd_uniform_latlon",
    "upwind_1st",
    "covariant_laplacian_cubed_sphere",
    "nn_diffusion_mpas",
    "ppm_reconstruction",
    "weno5_advection",
    "weno5_advection_2d",
    "flux_limiter_minmod",
    "flux_limiter_superbee",
    "divergence_arakawa_c",
    "transport_2d",
)

ALL_GRID_FAMILIES = (
    "cartesian",
    "latlon",
    "cubed_sphere",
    "mpas",
    "duo",
    "vertical",
    "arakawa",
)

EARTH_R = 1.0  # unit sphere is fine for visualizations


# ---------------------------------------------------------------------------
# Grid family visualizations
# ---------------------------------------------------------------------------


def plot_cartesian(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    nx, ny = 8, 6
    xs = np.linspace(0, 1, nx + 1)
    ys = np.linspace(0, 0.75, ny + 1)
    for x in xs:
        ax.plot([x, x], [ys[0], ys[-1]], color="#1f4f3a", lw=0.7)
    for y in ys:
        ax.plot([xs[0], xs[-1]], [y, y], color="#1f4f3a", lw=0.7)
    cx = 0.5 * (xs[:-1] + xs[1:])
    cy = 0.5 * (ys[:-1] + ys[1:])
    CX, CY = np.meshgrid(cx, cy)
    ax.scatter(CX, CY, color="#b58a00", s=18, zorder=5, label="cell centers")
    ax.set_aspect("equal")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 0.80)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Cartesian — 2D uniform rectilinear mesh (8×6)")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_latlon(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    nlon, nlat = 24, 12
    lons = np.linspace(-180, 180, nlon + 1)
    lats = np.linspace(-90, 90, nlat + 1)
    for lon in lons:
        ax.plot([lon, lon], [lats[0], lats[-1]], color="#1f4f3a", lw=0.6)
    for lat in lats:
        ax.plot([lons[0], lons[-1]], [lat, lat], color="#1f4f3a", lw=0.6)
    ax.set_aspect("equal")
    ax.set_xlim(-185, 185)
    ax.set_ylim(-95, 95)
    ax.set_xlabel("longitude (deg)")
    ax.set_ylabel("latitude (deg)")
    ax.set_title("Lat-Lon — regular 24×12 mesh (cell area ∝ cos φ)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_cubed_sphere(out: Path) -> None:
    """3D sphere with the 6 cubed-sphere panel boundaries overlaid."""
    fig = plt.figure(figsize=(5.5, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    # Sphere wireframe.
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(sx, sy, sz, color="#eef3f0", alpha=0.35, linewidth=0)

    # Six panel face-grid lines via gnomonic projection of a cube face.
    Nc = 6
    grid = np.linspace(-1, 1, Nc + 1)
    panels = []
    # +X, -X, +Y, -Y, +Z, -Z
    for axis in range(3):
        for sign in (+1, -1):
            face_pts = []
            for i in grid:
                row = []
                for j in grid:
                    p = [0.0, 0.0, 0.0]
                    p[axis] = sign
                    p[(axis + 1) % 3] = i
                    p[(axis + 2) % 3] = j
                    n = np.linalg.norm(p)
                    row.append([p[0] / n, p[1] / n, p[2] / n])
                face_pts.append(row)
            panels.append(np.array(face_pts))

    colors = ["#1f4f3a", "#b58a00", "#15315d", "#6c1d1d", "#5b3d09", "#1c4a31"]
    for face, col in zip(panels, colors):
        # iso-i lines
        for i in range(Nc + 1):
            ax.plot(face[i, :, 0], face[i, :, 1], face[i, :, 2], color=col, lw=0.6)
        # iso-j lines
        for j in range(Nc + 1):
            ax.plot(face[:, j, 0], face[:, j, 1], face[:, j, 2], color=col, lw=0.6)

    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    ax.set_title("Cubed-sphere — gnomonic C6 (6 panels × 6×6 cells)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_mpas(out: Path) -> None:
    """Quasi-uniform Voronoi cells on a sphere (illustrative; actual MPAS
    meshes have icosahedral structure with 12 pentagons + hexagons)."""
    from scipy.spatial import SphericalVoronoi
    rng = np.random.default_rng(7)
    n = 162
    pts = rng.normal(size=(n, 3))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    sv = SphericalVoronoi(pts, radius=1.0, center=np.array([0, 0, 0]))
    sv.sort_vertices_of_regions()

    fig = plt.figure(figsize=(5.5, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(sx, sy, sz, color="#f4f8f6", alpha=0.30, linewidth=0)

    for region in sv.regions:
        verts = sv.vertices[region]
        verts = np.vstack([verts, verts[0]])
        ax.plot(verts[:, 0], verts[:, 1], verts[:, 2], color="#1f4f3a", lw=0.5)
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], color="#b58a00", s=4)

    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    ax.set_title("MPAS — quasi-uniform spherical Voronoi mesh (162 cells)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def _icosahedron_vertices() -> np.ndarray:
    phi = (1 + 5 ** 0.5) / 2
    v = np.array(
        [
            [-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
            [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
            [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1],
        ],
        dtype=float,
    )
    v /= np.linalg.norm(v, axis=1)[:, None]
    return v


def _icosahedron_faces() -> list[tuple[int, int, int]]:
    return [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]


def plot_duo(out: Path) -> None:
    """Triangular icosahedral mesh (level 1: subdivide each face once)."""
    base = _icosahedron_vertices()
    base_faces = _icosahedron_faces()

    # Subdivide each triangle once to ~80 faces.
    verts = list(base)
    faces: list[tuple[int, int, int]] = []
    midpoint_cache: dict[tuple[int, int], int] = {}

    def midpoint(i: int, j: int) -> int:
        key = (min(i, j), max(i, j))
        if key in midpoint_cache:
            return midpoint_cache[key]
        m = (verts[i] + verts[j]) / 2
        m /= np.linalg.norm(m)
        verts.append(m)
        idx = len(verts) - 1
        midpoint_cache[key] = idx
        return idx

    for a, b, c in base_faces:
        ab = midpoint(a, b)
        bc = midpoint(b, c)
        ca = midpoint(c, a)
        faces.extend([(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)])

    V = np.array(verts)
    fig = plt.figure(figsize=(5.5, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(sx, sy, sz, color="#f4f8f6", alpha=0.30, linewidth=0)
    for tri in faces:
        loop = list(tri) + [tri[0]]
        pts = V[loop]
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color="#1f4f3a", lw=0.5)
    ax.scatter(V[:, 0], V[:, 1], V[:, 2], color="#b58a00", s=6)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    ax.set_title("Duo — geodesic triangular mesh (icosahedral, level 1)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_vertical(out: Path) -> None:
    """Stacked vertical levels (sigma-like spacing)."""
    fig, ax = plt.subplots(figsize=(4.0, 5.5))
    n = 12
    # Sigma-like: dense near the surface, sparse aloft.
    sigma = 1 - np.linspace(0, 1, n + 1) ** 1.6
    sigma = sigma[::-1]  # 0=top, 1=surface
    z_top = 20.0
    z_edges = (1 - sigma) * z_top
    for z in z_edges:
        ax.plot([0, 1], [z, z], color="#1f4f3a", lw=0.8)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    ax.scatter(np.full_like(z_centers, 0.5), z_centers, color="#b58a00", s=24, zorder=5)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.5, z_top + 0.5)
    ax.set_xticks([])
    ax.set_xlabel("(horizontal)")
    ax.set_ylabel("z (km)")
    ax.set_title("Vertical — 12 levels, stretched")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_arakawa(out: Path) -> None:
    """Arakawa A/B/C/D-grid stagger comparison (4 mini-cells)."""
    fig, axs = plt.subplots(1, 4, figsize=(11.5, 3.4))
    cell_pts = [(0, 0), (1, 0), (0, 1), (1, 1)]
    for ax, name in zip(axs, ("A", "B", "C", "D")):
        for x in (0, 1):
            ax.plot([x, x], [0, 1], color="#1f4f3a", lw=0.8)
        for y in (0, 1):
            ax.plot([0, 1], [y, y], color="#1f4f3a", lw=0.8)
        # Cell center mass marker
        ax.scatter([0.5], [0.5], color="#b58a00", s=80, marker="s", label="h (mass)")
        if name == "A":
            ax.scatter([0.5], [0.5], color="#15315d", s=60, marker=">", label="u")
            ax.scatter([0.5], [0.5], color="#6c1d1d", s=60, marker="^", label="v")
        elif name == "B":
            ax.scatter([1, 1, 0, 0], [0, 1, 0, 1], color="#15315d", s=60, marker=">",
                       label="u", zorder=5)
            ax.scatter([1, 1, 0, 0], [0, 1, 0, 1], color="#6c1d1d", s=60, marker="^",
                       label="v", zorder=6)
        elif name == "C":
            ax.scatter([0, 1], [0.5, 0.5], color="#15315d", s=70, marker=">", label="u (face-x)")
            ax.scatter([0.5, 0.5], [0, 1], color="#6c1d1d", s=70, marker="^", label="v (face-y)")
        elif name == "D":
            ax.scatter([0.5, 0.5], [0, 1], color="#15315d", s=70, marker=">", label="u (face-y)")
            ax.scatter([0, 1], [0.5, 0.5], color="#6c1d1d", s=70, marker="^", label="v (face-x)")
        ax.set_title(f"Arakawa {name}-grid")
        ax.set_xlim(-0.15, 1.15)
        ax.set_ylim(-0.15, 1.15)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), fontsize=7,
                  ncol=2, frameon=False)
    fig.suptitle("Arakawa staggers — h, u, v locations on a single cell", y=0.99)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


GRID_PLOTTERS = {
    "cartesian": plot_cartesian,
    "latlon": plot_latlon,
    "cubed_sphere": plot_cubed_sphere,
    "mpas": plot_mpas,
    "duo": plot_duo,
    "vertical": plot_vertical,
    "arakawa": plot_arakawa,
}


# ---------------------------------------------------------------------------
# Stencil / coefficient diagrams
# ---------------------------------------------------------------------------


def _stencil_1d(out: Path, title: str, points: list[tuple[int, str]],
                axis_label: str = "$x$ (offset)") -> None:
    """1D stencil with annotated coefficients above each grid point."""
    fig, ax = plt.subplots(figsize=(7.5, 3.0))
    offsets = [p[0] for p in points]
    coeffs = [p[1] for p in points]
    xmin, xmax = min(offsets) - 1, max(offsets) + 1
    ax.axhline(0, color="#1f4f3a", lw=0.8, zorder=1)
    for x in range(xmin, xmax + 1):
        ax.scatter([x], [0], color="#cccccc", s=80, zorder=2)
        ax.text(x, -0.20, f"i{x:+d}" if x != 0 else "i", ha="center",
                va="top", fontsize=9, color="#666")
    for o, c in zip(offsets, coeffs):
        ax.scatter([o], [0], color="#b58a00", s=160, zorder=3,
                   edgecolor="#5b3d09")
        ax.annotate(c, xy=(o, 0), xytext=(o, 0.30), ha="center",
                    fontsize=11, color="#1f4f3a")
    ax.set_xlim(xmin - 0.5, xmax + 0.5)
    ax.set_ylim(-0.7, 0.9)
    ax.set_yticks([])
    ax.set_xticks(list(range(xmin, xmax + 1)))
    ax.set_xticklabels([])
    ax.set_xlabel(axis_label)
    ax.set_title(title)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def _stencil_2d(out: Path, title: str,
                points: list[tuple[int, int, str]],
                xlabel: str = "$\\xi$ offset",
                ylabel: str = "$\\eta$ offset") -> None:
    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs) - 1, max(xs) + 1
    ymin, ymax = min(ys) - 1, max(ys) + 1
    for i in range(xmin, xmax + 1):
        ax.axvline(i, color="#e8e8e8", lw=0.5, zorder=0)
    for j in range(ymin, ymax + 1):
        ax.axhline(j, color="#e8e8e8", lw=0.5, zorder=0)
    for i in range(xmin, xmax + 1):
        for j in range(ymin, ymax + 1):
            ax.scatter([i], [j], color="#cccccc", s=40, zorder=1)
    for x, y, label in points:
        ax.scatter([x], [y], color="#b58a00", s=180, zorder=3,
                   edgecolor="#5b3d09")
        ax.text(x, y + 0.30, label, ha="center", fontsize=9, color="#1f4f3a")
    ax.set_aspect("equal")
    ax.set_xlim(xmin - 0.5, xmax + 0.5)
    ax.set_ylim(ymin - 0.5, ymax + 0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(list(range(xmin, xmax + 1)))
    ax.set_yticks(list(range(ymin, ymax + 1)))
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_centered_2nd(out: Path) -> None:
    _stencil_1d(out,
        "centered_2nd_uniform — ∂u/∂x stencil",
        [(-1, "−1/(2 dx)"), (1, "+1/(2 dx)")],
    )


def stencil_centered_2nd_vertical(out: Path) -> None:
    _stencil_1d(out,
        "centered_2nd_uniform_vertical — ∂u/∂k stencil",
        [(-1, "−1/(2 h)"), (1, "+1/(2 h)")],
        axis_label="$k$ (vertical level offset)",
    )


def stencil_centered_2nd_latlon(out: Path) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.0))
    for ax, axis_name, points in (
        (axs[0], "lon", [(-1, "−1/(2 R cosφ dλ)"), (1, "+1/(2 R cosφ dλ)")]),
        (axs[1], "lat", [(-1, "−1/(2 R dφ)"), (1, "+1/(2 R dφ)")]),
    ):
        offsets = [p[0] for p in points]
        ax.axhline(0, color="#1f4f3a", lw=0.8, zorder=1)
        for x in (-2, -1, 0, 1, 2):
            ax.scatter([x], [0], color="#cccccc", s=70, zorder=2)
            ax.text(x, -0.20, f"i{x:+d}" if x != 0 else "i", ha="center",
                    va="top", fontsize=9, color="#666")
        for o, c in points:
            ax.scatter([o], [0], color="#b58a00", s=140, zorder=3,
                       edgecolor="#5b3d09")
            ax.annotate(c, xy=(o, 0), xytext=(o, 0.30), ha="center",
                        fontsize=10, color="#1f4f3a")
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-0.7, 0.9)
        ax.set_yticks([])
        ax.set_xticks([-2, -1, 0, 1, 2])
        ax.set_xticklabels([])
        ax.set_xlabel(f"{axis_name} offset")
        ax.set_title(f"axis = {axis_name}")
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    fig.suptitle("centered_2nd_uniform_latlon — ∂u/∂x stencil with metric corrections")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_upwind_1st(out: Path) -> None:
    _stencil_1d(out,
        "upwind_1st — ∂u/∂x stencil (u > 0 branch)",
        [(-1, "−1/dx"), (0, "+1/dx")],
    )


def stencil_transport_2d(out: Path) -> None:
    pts = [
        (0, 0, "−Σ_face inflow / A"),
        (1, 0, "−((c_E−|c_E|)/2)·dx_E / A"),
        (-1, 0, "+((c_W+|c_W|)/2)·dx_W / A"),
        (0, 1, "−((c_N−|c_N|)/2)·dy_N / A"),
        (0, -1, "+((c_S+|c_S|)/2)·dy_S / A"),
    ]
    _stencil_2d(out, "transport_2d — 5-point in-panel stencil (C-grid Courant)", pts)


def stencil_covariant_laplacian(out: Path) -> None:
    pts = [
        (0, 0, "−2(g^{ξξ}+g^{ηη})/h²"),
        (1, 0, "g^{ξξ}/h² + ∂(Jg^{ξξ})/∂ξ +\n∂(Jg^{ξη})/∂η"),
        (-1, 0, "g^{ξξ}/h² − …"),
        (0, 1, "g^{ηη}/h² + …"),
        (0, -1, "g^{ηη}/h² − …"),
        (1, 1, "+g^{ξη}/(2 h²)"),
        (-1, -1, "+g^{ξη}/(2 h²)"),
        (1, -1, "−g^{ξη}/(2 h²)"),
        (-1, 1, "−g^{ξη}/(2 h²)"),
    ]
    _stencil_2d(out, "covariant_laplacian_cubed_sphere — 9-point in-panel stencil", pts)


def stencil_nn_diffusion_mpas(out: Path) -> None:
    """Voronoi neighbor diagram (illustrative)."""
    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    angles = np.linspace(0, 2 * np.pi, 7)[:-1]
    nbr = np.column_stack([np.cos(angles), np.sin(angles)])
    # central + 6 hex neighbors
    ax.scatter([0], [0], color="#b58a00", s=240, zorder=4,
               edgecolor="#5b3d09", label="target cell c")
    for i, (x, y) in enumerate(nbr):
        ax.scatter([x], [y], color="#dec98e", s=160, zorder=3,
                   edgecolor="#5b3d09")
        ax.plot([0, x], [0, y], color="#1f4f3a", lw=0.8, zorder=2)
        ax.text(x * 1.18, y * 1.18, f"c$_{{{i+1}}}$", ha="center", va="center",
                fontsize=10)
    # Hex outline (Voronoi cell)
    poly = np.array(
        [(np.cos(a), np.sin(a)) for a in np.linspace(np.pi/6, 2*np.pi+np.pi/6, 7)]
    ) * 0.55
    ax.plot(poly[:, 0], poly[:, 1], color="#1f4f3a", lw=1.0)
    ax.set_aspect("equal")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("nn_diffusion_mpas — neighbor stencil over edges_on_cell\n"
                 "coeff(c→c$_k$) = dv_edge / (dc_edge · area_cell)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_ppm_reconstruction(out: Path) -> None:
    pts = [(-2, "support"), (-1, "+7/12"), (0, "i"), (1, "+7/12"), (2, "support")]
    fig, ax = plt.subplots(figsize=(8.0, 3.0))
    ax.axhline(0, color="#1f4f3a", lw=0.8)
    for x in range(-3, 4):
        ax.scatter([x], [0], color="#cccccc", s=70, zorder=2)
        ax.text(x, -0.20, f"i{x:+d}" if x != 0 else "i", ha="center",
                va="top", fontsize=9, color="#666")
    edge_pts = {-1: "−1/12", 0: "+7/12", 1: "+7/12", 2: "−1/12"}
    for o, c in edge_pts.items():
        ax.scatter([o], [0], color="#b58a00", s=160, zorder=3,
                   edgecolor="#5b3d09")
        ax.annotate(c, xy=(o, 0), xytext=(o, 0.30), ha="center",
                    fontsize=10, color="#1f4f3a")
    # mark edge between i and i+1
    ax.annotate("q$_{i+1/2}$", xy=(0.5, 0), xytext=(0.5, -0.55),
                ha="center", fontsize=11, color="#6c1d1d",
                arrowprops=dict(arrowstyle="->", color="#6c1d1d"))
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-0.85, 0.9)
    ax.set_yticks([])
    ax.set_xticks(list(range(-3, 4)))
    ax.set_xticklabels([])
    ax.set_xlabel("$x$ (cell offset)")
    ax.set_title("ppm_reconstruction — 4th-order edge value (CW84 eq. 1.6)")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_weno5(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    sub_offsets = {
        "p₀": [-2, -1, 0],
        "p₁": [-1, 0, 1],
        "p₂": [0, 1, 2],
    }
    sub_y = {"p₀": 0.6, "p₁": 0.0, "p₂": -0.6}
    sub_color = {"p₀": "#15315d", "p₁": "#b58a00", "p₂": "#6c1d1d"}
    for x in range(-3, 4):
        for y in (-0.6, 0.0, 0.6):
            ax.scatter([x], [y], color="#eeeeee", s=40, zorder=1)
    for name, offs in sub_offsets.items():
        y = sub_y[name]
        for o in offs:
            ax.scatter([o], [y], color=sub_color[name], s=140, zorder=3,
                       edgecolor="#222")
        ax.text(3.6, y, f"{name}  (d = {('1/10','6/10','3/10')[list(sub_offsets).index(name)]})",
                fontsize=10, color=sub_color[name], va="center")
    ax.text(0.5, 1.10, "edge q$_{i+1/2}^L$",
            ha="center", color="#6c1d1d", fontsize=11)
    ax.axvline(0.5, color="#6c1d1d", lw=0.6, ls="--")
    for x in range(-3, 4):
        ax.text(x, -1.1, f"i{x:+d}" if x != 0 else "i", ha="center",
                va="top", fontsize=9, color="#666")
    ax.set_xlim(-3.5, 5.5)
    ax.set_ylim(-1.4, 1.3)
    ax.set_yticks([])
    ax.set_xticks(list(range(-3, 4)))
    ax.set_xticklabels([])
    ax.set_title("weno5_advection — three sub-stencils (left-biased branch)")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_weno5_2d(out: Path) -> None:
    """Cross-shaped 2D WENO5 stencil: dimension-by-dimension splitting reuses
    the 1D Jiang-Shu (1996) sub-stencils along x and y independently. Nine
    points along each axis cross at (i, j); the two arms together span 17 of
    the 25 cells in the 5×5 neighborhood."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for i in range(-3, 4):
        for j in range(-3, 4):
            ax.scatter([i], [j], color="#eeeeee", s=40, zorder=1)
    x_arm = [(-2, 0), (-1, 0), (0, 0), (1, 0), (2, 0)]
    y_arm = [(0, -2), (0, -1), (0, 1), (0, 2)]
    for x, y in x_arm:
        ax.scatter([x], [y], color="#15315d", s=170, zorder=3,
                   edgecolor="#0a1530")
    for x, y in y_arm:
        ax.scatter([x], [y], color="#6c1d1d", s=170, zorder=3,
                   edgecolor="#3a0a0a")
    ax.scatter([0], [0], color="#b58a00", s=180, zorder=4,
               edgecolor="#5b3d09")
    ax.axvline(0.5, color="#15315d", lw=0.6, ls="--",
               label="x face $x_{i+1/2}$")
    ax.axhline(0.5, color="#6c1d1d", lw=0.6, ls="--",
               label="y face $y_{j+1/2}$")
    for x in range(-3, 4):
        ax.text(x, -3.65, f"i{x:+d}" if x != 0 else "i", ha="center",
                fontsize=8, color="#666")
    for y in range(-3, 4):
        ax.text(-3.6, y, f"j{y:+d}" if y != 0 else "j", va="center",
                fontsize=8, color="#666")
    ax.set_xlim(-3.7, 3.7)
    ax.set_ylim(-3.7, 3.7)
    ax.set_aspect("equal")
    ax.set_xticks(list(range(-3, 4)))
    ax.set_yticks(list(range(-3, 4)))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title("weno5_advection_2d — dimension-by-dimension cross stencil")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def _limiter_curve(out: Path, title: str, phi):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    r = np.linspace(-0.5, 4, 600)
    ax.plot(r, [phi(ri) for ri in r], color="#1f4f3a", lw=2.0, label="φ(r)")
    # Sweby second-order TVD region (between phi=r and phi=1, etc.).
    rr = np.linspace(0, 4, 400)
    upper = np.minimum(2 * rr, np.minimum(2.0, np.maximum(rr, 1.0)))
    upper2 = np.where(rr < 1, 2 * rr, 2.0)
    lower = np.where(rr < 1, rr, 1.0)
    ax.fill_between(rr, lower, upper2, color="#dde7e0", alpha=0.5,
                    label="Sweby 2nd-order TVD region")
    ax.plot([0, 4], [0, 0], color="#888", lw=0.5)
    ax.plot([0, 0], [-0.3, 2.3], color="#888", lw=0.5)
    ax.scatter([1], [1], color="#b58a00", zorder=5)
    ax.annotate("φ(1) = 1", xy=(1, 1), xytext=(1.4, 1.4),
                arrowprops=dict(arrowstyle="->", color="#b58a00"),
                color="#5b3d09")
    ax.set_xlim(-0.5, 4)
    ax.set_ylim(-0.3, 2.3)
    ax.set_xlabel("slope ratio  r = (qᵢ − qᵢ₋₁) / (qᵢ₊₁ − qᵢ)")
    ax.set_ylabel("limiter  φ(r)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_limiter_minmod(out: Path) -> None:
    _limiter_curve(out,
        "flux_limiter_minmod — φ(r) = max(0, min(r, 1))",
        lambda r: max(0.0, min(r, 1.0)),
    )


def stencil_limiter_superbee(out: Path) -> None:
    _limiter_curve(out,
        "flux_limiter_superbee — φ(r) = max(0, min(2r,1), min(r,2))",
        lambda r: max(0.0, min(2 * r, 1.0), min(r, 2.0)),
    )


def interface_limiter_minmod(out: Path) -> None:
    """Slope-ratio interface stencil: q_{i-1}, q_i, q_{i+1} feeding face i+1/2."""
    fig, ax = plt.subplots(figsize=(7.8, 3.4))
    centers_x = [-1, 0, 1]
    center_labels = ["qᵢ₋₁", "qᵢ", "qᵢ₊₁"]
    face_x = 0.5
    # Cell baseline
    ax.axhline(0, color="#1f4f3a", lw=0.8, zorder=1)
    # Cell boundaries
    for x in (-1.5, -0.5, 0.5, 1.5):
        ax.plot([x, x], [-0.18, 0.18], color="#888", lw=0.6)
    # Cell centers
    for x, lbl in zip(centers_x, center_labels):
        ax.scatter([x], [0], color="#b58a00", s=160, zorder=3,
                   edgecolor="#5b3d09")
        ax.text(x, -0.32, lbl, ha="center", va="top", fontsize=10,
                color="#1f4f3a")
    # Highlight the i+1/2 face
    ax.scatter([face_x], [0], color="#6c1d1d", s=180, marker="D",
               zorder=4, edgecolor="#3a0a0a", label="face i+1/2")
    ax.annotate("Fᵢ₊₁/₂  (limited flux emitted here)",
                xy=(face_x, 0), xytext=(face_x, 0.55),
                ha="center", fontsize=10, color="#6c1d1d",
                arrowprops=dict(arrowstyle="->", color="#6c1d1d"))
    # Slope-ratio annotation
    ax.annotate("", xy=(0, 0.22), xytext=(-1, 0.22),
                arrowprops=dict(arrowstyle="<->", color="#15315d", lw=0.9))
    ax.text(-0.5, 0.30, "qᵢ − qᵢ₋₁  (upwind slope)", ha="center",
            fontsize=9, color="#15315d")
    ax.annotate("", xy=(1, 0.22), xytext=(0, 0.22),
                arrowprops=dict(arrowstyle="<->", color="#5b3d09", lw=0.9))
    ax.text(0.5, 0.30, "qᵢ₊₁ − qᵢ  (downwind slope)", ha="center",
            fontsize=9, color="#5b3d09")
    ax.text(0.0, -0.78,
            "r = (qᵢ − qᵢ₋₁) / (qᵢ₊₁ − qᵢ)     "
            "→     φ(r) = max(0, min(r, 1))     "
            "→     Fᵢ₊₁/₂ = F_low + φ(r) · (F_high − F_low)",
            ha="center", fontsize=9.5, color="#1f4f3a")
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.0, 0.85)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title("flux_limiter_minmod — slope-ratio interface stencil "
                 "(positive-velocity branch)")
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def stencil_arakawa_divergence(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    # Single C-grid cell. cell-center at (0.5, 0.5).
    for x in (0, 1):
        ax.plot([x, x], [0, 1], color="#1f4f3a", lw=0.8)
    for y in (0, 1):
        ax.plot([0, 1], [y, y], color="#1f4f3a", lw=0.8)
    ax.scatter([0.5], [0.5], color="#b58a00", s=180, marker="s",
               edgecolor="#5b3d09", zorder=5, label="div emitted (cell_center)")
    ax.scatter([0, 1], [0.5, 0.5], color="#15315d", s=140, marker=">",
               edgecolor="#0a1530", zorder=5, label="F$_x$ (face_x)")
    ax.scatter([0.5, 0.5], [0, 1], color="#6c1d1d", s=140, marker="^",
               edgecolor="#3a0a0a", zorder=5, label="F$_y$ (face_y)")
    ax.annotate("−1/dx", xy=(0, 0.5), xytext=(-0.25, 0.55), fontsize=10,
                color="#15315d")
    ax.annotate("+1/dx", xy=(1, 0.5), xytext=(1.05, 0.55), fontsize=10,
                color="#15315d")
    ax.annotate("−1/dy", xy=(0.5, 0), xytext=(0.55, -0.15), fontsize=10,
                color="#6c1d1d")
    ax.annotate("+1/dy", xy=(0.5, 1), xytext=(0.55, 1.07), fontsize=10,
                color="#6c1d1d")
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylim(-0.3, 1.3)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("divergence_arakawa_c — C-grid 4-point stencil")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3,
              fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


RULE_STENCIL_PLOTTERS = {
    "centered_2nd_uniform": stencil_centered_2nd,
    "centered_2nd_uniform_vertical": stencil_centered_2nd_vertical,
    "centered_2nd_uniform_latlon": stencil_centered_2nd_latlon,
    "upwind_1st": stencil_upwind_1st,
    "covariant_laplacian_cubed_sphere": stencil_covariant_laplacian,
    "nn_diffusion_mpas": stencil_nn_diffusion_mpas,
    "ppm_reconstruction": stencil_ppm_reconstruction,
    "weno5_advection": stencil_weno5,
    "weno5_advection_2d": stencil_weno5_2d,
    "flux_limiter_minmod": stencil_limiter_minmod,
    "flux_limiter_superbee": stencil_limiter_superbee,
    "divergence_arakawa_c": stencil_arakawa_divergence,
    "transport_2d": stencil_transport_2d,
}


# Optional secondary diagrams keyed by suffix. A rule appears here only if it
# benefits from a second image (e.g. limiters: φ(r) curve + grid-space
# slope-ratio stencil). File naming is rules/<rule>-<suffix>.png.
RULE_EXTRA_PLOTTERS: dict[str, dict[str, "callable"]] = {
    "flux_limiter_minmod": {"interface": interface_limiter_minmod},
}


# ---------------------------------------------------------------------------
# Convergence plots (active rules only)
#
# These re-implement the rule + manufactured solution from the convergence
# fixture under ``discretizations/<family>/<rule>/fixtures/convergence/`` so
# the doc page can show the empirical order of convergence. The walker
# remains the authoritative oracle — these are illustrative.
# ---------------------------------------------------------------------------


def _convergence_plot(out: Path, title: str, ns, errs, expected_order: float) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(ns, errs, "o-", color="#b58a00", lw=1.5, markersize=7,
              markeredgecolor="#5b3d09", label="observed L∞ error")
    ref = errs[0] * (ns[0] / np.array(ns, dtype=float)) ** expected_order
    ax.loglog(ns, ref, "--", color="#1f4f3a", lw=1.2,
              label=f"reference slope −{expected_order:.0f}")
    # Estimate observed order from fit.
    p = np.polyfit(np.log(ns), np.log(errs), 1)
    ax.set_xlabel("number of cells N")
    ax.set_ylabel("L∞ error")
    ax.set_title(title + f"\nfit slope = {p[0]:+.2f}")
    ax.grid(True, which="both", color="#eee", lw=0.5)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def convergence_centered_2nd(out: Path) -> None:
    """sin(2πx) on [0,1] periodic, cell-centered samples,
    centered_2nd_uniform applied to grad."""
    ns = [16, 32, 64, 128, 256]
    errs = []
    for n in ns:
        dx = 1.0 / n
        x = (np.arange(n) + 0.5) * dx
        u = np.sin(2 * np.pi * x)
        # Periodic centered diff on cell-centered values
        du_num = (np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)
        du_exact = 2 * np.pi * np.cos(2 * np.pi * x)
        errs.append(np.max(np.abs(du_num - du_exact)))
    _convergence_plot(out, "centered_2nd_uniform — empirical convergence",
                      ns, errs, expected_order=2.0)


def convergence_centered_2nd_vertical(out: Path) -> None:
    """Same kernel as horizontal — vertical only relabels axis."""
    ns = [16, 32, 64, 128, 256]
    errs = []
    for n in ns:
        dz = 1.0 / n
        z = (np.arange(n) + 0.5) * dz
        u = np.cos(2 * np.pi * z)
        du_num = (np.roll(u, -1) - np.roll(u, 1)) / (2 * dz)
        du_exact = -2 * np.pi * np.sin(2 * np.pi * z)
        errs.append(np.max(np.abs(du_num - du_exact)))
    _convergence_plot(out,
                      "centered_2nd_uniform_vertical — empirical convergence",
                      ns, errs, expected_order=2.0)


def convergence_upwind_1st(out: Path) -> None:
    """Upwind 1st on grad. u > 0 branch.  manufactured: sin(2πx)."""
    ns = [16, 32, 64, 128, 256, 512]
    errs = []
    for n in ns:
        dx = 1.0 / n
        x = (np.arange(n) + 0.5) * dx
        u = np.sin(2 * np.pi * x)
        # backward (u > 0): du[i] = (u[i] - u[i-1]) / dx
        du_num = (u - np.roll(u, 1)) / dx
        du_exact = 2 * np.pi * np.cos(2 * np.pi * x)
        errs.append(np.max(np.abs(du_num - du_exact)))
    _convergence_plot(out, "upwind_1st — empirical convergence",
                      ns, errs, expected_order=1.0)


def convergence_ppm_reconstruction(out: Path) -> None:
    """Unlimited PPM (CW84 §1, eqs. 1.6, 1.7, 1.10) on sin(2πx) on [0,1]
    periodic, cell-averaged inputs, evaluated at subcell ξ ∈ {0.1, …, 0.9}.
    Mirrors the fixture under
    ``discretizations/finite_volume/ppm_reconstruction/fixtures/convergence``.
    """
    ns = [16, 32, 64, 128, 256]
    xis = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    errs = []
    for n in ns:
        dx = 1.0 / n
        # Cell centers and exact cell averages of sin(2πx) on [x-dx/2, x+dx/2].
        x = (np.arange(n) + 0.5) * dx
        q = (np.cos(2 * np.pi * (x - 0.5 * dx))
             - np.cos(2 * np.pi * (x + 0.5 * dx))) / (2 * np.pi * dx)
        # 4th-order edge value at i+1/2: (-q_{i-1} + 7 q_i + 7 q_{i+1} - q_{i+2}) / 12.
        q_edge = (-np.roll(q, 1) + 7 * q + 7 * np.roll(q, -1)
                  - np.roll(q, -2)) / 12.0
        a_L = np.roll(q_edge, 1)   # q_{i-1/2}
        a_R = q_edge               # q_{i+1/2}
        da = a_R - a_L
        a6 = 6.0 * (q - 0.5 * (a_L + a_R))
        # Evaluate parabola at each ξ for each cell, compare to exact sin.
        max_err = 0.0
        for xi in xis:
            recon = a_L + xi * (da + a6 * (1.0 - xi))
            exact = np.sin(2 * np.pi * (x - 0.5 * dx + xi * dx))
            max_err = max(max_err, float(np.max(np.abs(recon - exact))))
        errs.append(max_err)
    _convergence_plot(out, "ppm_reconstruction — empirical convergence",
                      ns, errs, expected_order=3.0)


def convergence_centered_2nd_latlon(out: Path) -> None:
    """Y_{2,0} spherical harmonic on the unit sphere, lon-independent so the
    lon stencil is exact and the test isolates the lat-axis order. Mirrors the
    Layer-B fixture under
    discretizations/finite_difference/centered_2nd_uniform_latlon/fixtures/convergence/.
    """
    R = EARTH_R
    ns = [16, 32, 64, 128]
    errs = []
    for n in ns:
        # Cell-centered lat grid on [-π/2, π/2] with constant dφ = π/n.
        dphi = math.pi / n
        phi = (np.arange(n) + 0.5) * dphi - math.pi / 2.0
        u = 3.0 * np.sin(phi) ** 2 - 1.0  # ∝ Y_{2,0}(φ)
        # Centered lat stencil; interior cells only (poles excluded).
        du_num = (u[2:] - u[:-2]) / (2.0 * R * dphi)
        du_exact = 3.0 * np.sin(2.0 * phi[1:-1]) / R
        errs.append(np.max(np.abs(du_num - du_exact)))
    _convergence_plot(out,
                      "centered_2nd_uniform_latlon — empirical convergence",
                      ns, errs, expected_order=2.0)


def convergence_weno5_advection_2d(out: Path) -> None:
    """Dimension-by-dimension Jiang-Shu (1996) WENO5 reconstruction on the
    phase-shifted 2D sine u(x,y) = sin(2π x + 1.0)·sin(2π y + 0.5) on
    [0,1]² periodic, cell-averaged inputs. Mirrors the Layer-C fixture under
    discretizations/finite_volume/weno5_advection_2d/fixtures/integration/
    smooth_2d_convergence.esm. Asymptotic order is sub-5th (~4.5–4.7) with
    ε = 1e-6 nonlinear-weight regularisation; threshold 4.4 across n ∈
    {16,32,64,128} per dsc-5od acceptance."""
    ns = [16, 32, 64, 128]
    phi_x, phi_y = 1.0, 0.5
    eps = 1.0e-6
    errs = []

    def cell_avg_1d(a, b, phi):
        return (math.cos(2 * math.pi * a + phi)
                - math.cos(2 * math.pi * b + phi)) / (2 * math.pi * (b - a))

    def weno5_left_face(q):
        n = len(q)
        out = np.empty(n)
        for i in range(n):
            qm2 = q[(i - 2) % n]; qm1 = q[(i - 1) % n]; q0 = q[i]
            qp1 = q[(i + 1) % n]; qp2 = q[(i + 2) % n]
            p0 = (1 / 3) * qm2 + (-7 / 6) * qm1 + (11 / 6) * q0
            p1 = (-1 / 6) * qm1 + (5 / 6) * q0 + (1 / 3) * qp1
            p2 = (1 / 3) * q0 + (5 / 6) * qp1 + (-1 / 6) * qp2
            b0 = (13 / 12) * (qm2 - 2 * qm1 + q0) ** 2 \
                + 0.25 * (qm2 - 4 * qm1 + 3 * q0) ** 2
            b1 = (13 / 12) * (qm1 - 2 * q0 + qp1) ** 2 \
                + 0.25 * (qm1 - qp1) ** 2
            b2 = (13 / 12) * (q0 - 2 * qp1 + qp2) ** 2 \
                + 0.25 * (3 * q0 - 4 * qp1 + qp2) ** 2
            a0 = 0.1 / (eps + b0) ** 2
            a1 = 0.6 / (eps + b1) ** 2
            a2 = 0.3 / (eps + b2) ** 2
            s = a0 + a1 + a2
            out[i] = (a0 / s) * p0 + (a1 / s) * p1 + (a2 / s) * p2
        return out

    for n in ns:
        h = 1.0 / n
        x_avg = np.array([cell_avg_1d(i * h, (i + 1) * h, phi_x)
                          for i in range(n)])
        y_avg = np.array([cell_avg_1d(j * h, (j + 1) * h, phi_y)
                          for j in range(n)])
        qbar = np.outer(x_avg, y_avg)  # qbar[i, j]
        x_face = np.array([math.sin(2 * math.pi * (i + 1) * h + phi_x)
                           for i in range(n)])
        y_face = np.array([math.sin(2 * math.pi * (j + 1) * h + phi_y)
                           for j in range(n)])
        err = 0.0
        for j in range(n):
            qhat = weno5_left_face(qbar[:, j])
            truth = x_face * y_avg[j]
            err = max(err, float(np.max(np.abs(qhat - truth))))
        for i in range(n):
            qhat = weno5_left_face(qbar[i, :])
            truth = x_avg[i] * y_face
            err = max(err, float(np.max(np.abs(qhat - truth))))
        errs.append(err)
    _convergence_plot(out,
                      "weno5_advection_2d — empirical convergence",
                      ns, errs, expected_order=5.0)


CONVERGENCE_PLOTTERS = {
    "centered_2nd_uniform": convergence_centered_2nd,
    "centered_2nd_uniform_vertical": convergence_centered_2nd_vertical,
    "centered_2nd_uniform_latlon": convergence_centered_2nd_latlon,
    "ppm_reconstruction": convergence_ppm_reconstruction,
    "upwind_1st": convergence_upwind_1st,
    "weno5_advection_2d": convergence_weno5_advection_2d,
}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def render_all(out_root: Path, what: str = "all") -> None:
    grids_dir = out_root / "grids"
    rules_dir = out_root / "rules"
    _ensure_dir(grids_dir)
    _ensure_dir(rules_dir)

    if what in ("all", "grids"):
        for fam in ALL_GRID_FAMILIES:
            target = grids_dir / f"{fam}.png"
            print(f"  grids/{fam}.png", file=sys.stderr)
            GRID_PLOTTERS[fam](target)

    if what in ("all", "rules"):
        for rule in ALL_RULES:
            target = rules_dir / f"{rule}-stencil.png"
            print(f"  rules/{rule}-stencil.png", file=sys.stderr)
            RULE_STENCIL_PLOTTERS[rule](target)
            for suffix, plotter in RULE_EXTRA_PLOTTERS.get(rule, {}).items():
                extra = rules_dir / f"{rule}-{suffix}.png"
                print(f"  rules/{rule}-{suffix}.png", file=sys.stderr)
                plotter(extra)
        for rule in sorted(APPLICABLE):
            target = rules_dir / f"{rule}-convergence.png"
            print(f"  rules/{rule}-convergence.png", file=sys.stderr)
            CONVERGENCE_PLOTTERS[rule](target)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="docs/static/plots",
                    help="Output root (default: docs/static/plots)")
    ap.add_argument("--what", choices=("all", "grids", "rules"), default="all")
    args = ap.parse_args()
    render_all(Path(args.out), what=args.what)
    return 0


if __name__ == "__main__":
    sys.exit(main())
