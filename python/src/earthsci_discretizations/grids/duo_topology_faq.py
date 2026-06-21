"""DUO primal triangular-mesh topology via the landed ESS value-invention engine
(``earthsci_toolkit`` relational primitives).

Python mirror of ``src/topology_faq.jl`` (D1a). Declarative companion:
``discretizations/grids/duo/rules/primal_topology.esm`` expresses this chain as a
value-invention FAQ document (RFC ``semiring-faq-unified-ir`` §7.3: enumerate edges
from the face→vertex relation, ``skolem`` on the sorted vertex-pair key,
``distinct`` dedup, ``rank`` dense ids, then the inversion joins for
``cell_neighbors`` and the group-by for ``vertex_faces``).

This module is the thin **evaluation bridge**: it shapes a triangular ``faces``
matrix into the relational inputs and routes every dedup / rank / join through
ESS's relational primitives. ESD hosts no shadow relational engine — the
determinism contract (and the cross-binding byte-identity of the output) lives
entirely in ``earthsci_toolkit`` (CONFORMANCE_SPEC §5.5; AGENTS.md single-pathway
invariant). The only ESD-side logic is *authoring* which vertex pairs are a
triangle's edges.

Edge numbering is the ESS canonical (sorted) order produced by ``distinct`` —
the cross-binding byte-stable numbering, replacing the imperative
``_build_connectivity`` dict-insertion order (binding-divergent). ``cell_neighbors``
(cell×slot) and ``vertex_faces`` are byte-identical regardless of edge numbering.
Dense edge ids use the 0-based Python convention (``rank`` default ``base=0``); the
absolute base is an internal join key and does not affect ``cell_neighbors`` (cell
indices) or the canonical edge order.
"""

from __future__ import annotations

import numpy as np
from earthsci_toolkit import distinct, equijoin, rank, skolem_edge

__all__ = ["primal_topology_faq"]


def primal_topology_faq(faces: np.ndarray, n_vertices: int) -> dict:
    """Materialize DUO primal unstructured connectivity from a ``(3, Nc)``
    triangular ``faces`` matrix (0-based vertex ids) by evaluating the
    value-invention FAQ through the landed ESS relational engine.

    Returns a dict with the four primal-topology arrays, matching the
    :class:`DuoGrid` field conventions:

    * ``edges`` ``(2, Ne)`` int64 — sorted vertex pairs ``(a < b)``, emitted in the
      canonical (sorted) order.
    * ``edges_on_face`` ``(3, Nc)`` int64 — 0-based dense id of edge ``k`` of face
      ``c``, where slot ``k`` is the edge OPPOSITE vertex ``k`` (``k=0→(v1,v2)``,
      ``k=1→(v2,v0)``, ``k=2→(v0,v1)``), the same slot convention as
      ``cell_neighbors``.
    * ``cell_neighbors`` ``(3, Nc)`` int64 — cell sharing edge ``k`` of face ``c``;
      ``-1`` for a boundary edge (none on a closed icosahedral mesh).
    * ``vertex_faces`` ``list[list[int]]`` — sorted face ids incident on each vertex.

    The dedup (``distinct``), dense renumbering (``rank``), and connectivity
    inversion (``equijoin``) are computed by ``earthsci_toolkit``, so the result is
    byte-identical to the Rust / Julia bindings' engine by ESS's cross-binding
    determinism contract.
    """
    if faces.shape[0] != 3:
        raise ValueError(
            f"primal_topology_faq: faces must be (3, Nc) triangular; got {faces.shape}"
        )
    if n_vertices < 0:
        raise ValueError(
            f"primal_topology_faq: n_vertices must be >= 0, got {n_vertices}"
        )
    Nc = int(faces.shape[1])

    # The three undirected edges of triangle ``c``, indexed by the local slot ``k``
    # OPPOSITE vertex ``k`` (k=0→(v1,v2), k=1→(v2,v0), k=2→(v0,v1)) — the slot
    # convention shared by ``cell_neighbors``/``edges_on_face`` and matched to the
    # imperative ``_build_connectivity``.
    def _face_edges(c: int):
        v0, v1, v2 = int(faces[0, c]), int(faces[1, c]), int(faces[2, c])
        return ((v1, v2), (v2, v0), (v0, v1))

    # --- (1) edges: skolem on the sorted vertex pair + distinct + rank ----------
    edge_keys = []
    for c in range(Nc):
        for a, b in _face_edges(c):
            edge_keys.append(skolem_edge(a, b))
    edge_order = distinct(edge_keys)  # sorted unique (min, max) pairs
    ranking = rank(edge_keys, base=0)  # 0-based dense ids (Python native)
    Ne = len(edge_order)
    edges = np.empty((2, Ne), dtype=np.int64)
    for e, ab in enumerate(edge_order):
        edges[0, e] = ab[0]
        edges[1, e] = ab[1]

    # --- (2) edges_on_face: dense edge id of each of a face's three edges -------
    edges_on_face = np.empty((3, Nc), dtype=np.int64)
    for c in range(Nc):
        for k, (a, b) in enumerate(_face_edges(c)):
            edges_on_face[k, c] = ranking[skolem_edge(a, b)]

    # --- (3) cell_neighbors: the OTHER cell across each edge --------------------
    # Inversion join (ESS equijoin — the connectivity-inversion primitive):
    # self-join the (dense edge id, cell, slot) incidence on the edge id; each pair
    # (c1, c2) with c1 != c2 on a shared edge is the neighbour across that slot.
    incidence = []  # (dense edge id, cell, slot k)
    for c in range(Nc):
        for k, (a, b) in enumerate(_face_edges(c)):
            incidence.append((ranking[skolem_edge(a, b)], c, k))
    cell_neighbors = np.full((3, Nc), -1, dtype=np.int64)
    joined = equijoin(incidence, incidence, on_left=lambda t: t[0], on_right=lambda t: t[0])
    for left, right in joined:
        c1, k1 = left[1], left[2]
        c2 = right[1]
        if c1 == c2:
            continue
        if cell_neighbors[k1, c1] != -1 and cell_neighbors[k1, c1] != c2:
            raise AssertionError(
                f"primal_topology_faq: non-manifold edge — cell {c1} slot {k1} "
                "shares an edge with >2 cells"
            )
        cell_neighbors[k1, c1] = c2

    # --- (4) vertex_faces: faces incident on each vertex, sorted ----------------
    # Group-by over the face→vertex relation: ``distinct`` over (vertex, face) rows
    # sorts by (vertex, face), so each vertex's faces land contiguous and ascending.
    vf_rows = []
    for c in range(Nc):
        for k in range(3):
            vf_rows.append((int(faces[k, c]), c))
    vf_sorted = distinct(vf_rows)
    vertex_faces: list[list[int]] = [[] for _ in range(int(n_vertices))]
    for v, f in vf_sorted:
        vertex_faces[v].append(f)

    return {
        "edges": edges,
        "edges_on_face": edges_on_face,
        "cell_neighbors": cell_neighbors,
        "vertex_faces": vertex_faces,
    }
