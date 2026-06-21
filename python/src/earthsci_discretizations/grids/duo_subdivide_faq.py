"""DUO icosahedral mesh construction with every vertex COORDINATE routed through
the ESS front-door evaluator (``eval_coeff``), preserving the EXACT
combinatorial ORDER of the imperative ``_subdivide_icosahedron`` (:mod:`duo`).

Python mirror of ``src/subdivide_faq.jl`` (D3, esd-ohd / W1). Declarative
companions: ``discretizations/grids/duo/faq/subdivide_fold.esm`` and
``subdivide_refine.esm``. The seed ``vert_raw / sqrt(seed_sq)`` normalize and the
midpoint ``(a+b) / sqrt(Σ (a+b)²)`` average are the geometry FAQs this file
evaluates through ``eval_coeff`` (``earthsci_discretizations.rules`` — the single
arithmetic pathway; no shadow evaluator, per AGENTS.md / GRIDS_API §4.3).

Front-door coordinates, imperative encounter-order numbering
------------------------------------------------------------
Only the **coordinates** go through the front door: every vertex value (the 12
seed vertices and every refined midpoint) is computed by ``eval_coeff`` on the
SAME ASTs the Julia binding drives. The **numbering**, however, is the imperative
encounter order — 12 base vertices in table order, then midpoints minted in
``dict``-cache ENCOUNTER order as faces are iterated — NOT the canonical
sorted/skolem/rank order serialized in ``subdivide_refine.esm``.

The production path must keep the imperative order because downstream
dual-conformance asserts the per-vertex Voronoi-dual arrays equal the imperative
``_duo_voronoi_dual`` EXACTLY, and the dual is indexed by primal vertex — so the
primal vertex numbering must match ``_subdivide_icosahedron`` column-for-column.
``eval_coeff`` on these ASTs is bit-identical to the imperative
``_icosahedron_vertices`` / ``midpoint`` arithmetic (squares are products, the
3-arg ``+`` folds left-associatively, the single divide-by-``sqrt`` is the same
op), so the output is byte-identical AND same-order as
``_subdivide_icosahedron(level)`` for ``float64``.
"""

from __future__ import annotations

import math

import numpy as np

from ..rules import eval_coeff

__all__ = ["duo_subdivide_faq"]


# Raw golden-ratio icosahedron vertices, pre-normalization. Transcribed to match
# ``_icosahedron_vertices``'s ``raw`` literal table (:mod:`duo`) row-for-row;
# phi = (1 + sqrt(5))/2.
_FAQ_PHI = (1.0 + math.sqrt(5.0)) / 2.0
_FAQ_VERT_RAW = (
    (0.0, 1.0, _FAQ_PHI),
    (0.0, -1.0, _FAQ_PHI),
    (0.0, 1.0, -_FAQ_PHI),
    (0.0, -1.0, -_FAQ_PHI),
    (1.0, _FAQ_PHI, 0.0),
    (-1.0, _FAQ_PHI, 0.0),
    (1.0, -_FAQ_PHI, 0.0),
    (-1.0, -_FAQ_PHI, 0.0),
    (_FAQ_PHI, 0.0, 1.0),
    (_FAQ_PHI, 0.0, -1.0),
    (-_FAQ_PHI, 0.0, 1.0),
    (-_FAQ_PHI, 0.0, -1.0),
)


# --- Seed normalize-to-sphere FAQ AST ----------------------------------------
# coord_d = r$d / sqrt(r1*r1 + r2*r2 + r3*r3).  Squares are products (``*``), NOT
# ``^``, and the 3-arg ``+`` folds left-associatively — the exact float-op
# sequence of the imperative ``n = sqrt(v0**2 + v1**2 + v2**2); v[d]/n``.
def _faq_rsq(name: str) -> dict:
    return {"op": "*", "args": [name, name]}


_FAQ_SEED_NORM = {
    "op": "sqrt",
    "args": [{"op": "+", "args": [_faq_rsq("r1"), _faq_rsq("r2"), _faq_rsq("r3")]}],
}
_FAQ_SEED_COMP = (
    {"op": "/", "args": ["r1", _FAQ_SEED_NORM]},
    {"op": "/", "args": ["r2", _FAQ_SEED_NORM]},
    {"op": "/", "args": ["r3", _FAQ_SEED_NORM]},
)


# --- Midpoint normalize-to-sphere FAQ AST ------------------------------------
# coord_d = (a$d + b$d) / sqrt(Σ_d (a$d + b$d)*(a$d + b$d)).  Squares are products
# so the float ops match the imperative ``midpoint`` bit-for-bit.
def _faq_sum(d: int) -> dict:
    return {"op": "+", "args": [f"a{d}", f"b{d}"]}


def _faq_sq(d: int) -> dict:
    return {"op": "*", "args": [_faq_sum(d), _faq_sum(d)]}


_FAQ_MID_NORM = {
    "op": "sqrt",
    "args": [{"op": "+", "args": [_faq_sq(1), _faq_sq(2), _faq_sq(3)]}],
}
_FAQ_MID_COMP = (
    {"op": "/", "args": [_faq_sum(1), _FAQ_MID_NORM]},
    {"op": "/", "args": [_faq_sum(2), _FAQ_MID_NORM]},
    {"op": "/", "args": [_faq_sum(3), _FAQ_MID_NORM]},
)


def _faq_seed_vertices() -> list[tuple[float, float, float]]:
    """The 12 base icosahedron vertices on the unit sphere, in ``_FAQ_VERT_RAW``
    table order (= ``_icosahedron_vertices`` order), each computed via the seed
    normalize FAQ (``eval_coeff``). Bit-identical to ``_icosahedron_vertices()``
    columns for ``float64``. Values are full-precision Python floats (``float64``);
    the final dtype cast happens when the vertex array is materialized, mirroring
    the imperative path which subdivides in ``float64`` regardless of output dtype.
    """
    verts: list[tuple[float, float, float]] = []
    for r1, r2, r3 in _FAQ_VERT_RAW:
        b = {"r1": float(r1), "r2": float(r2), "r3": float(r3)}
        verts.append(
            (
                eval_coeff(_FAQ_SEED_COMP[0], b),
                eval_coeff(_FAQ_SEED_COMP[1], b),
                eval_coeff(_FAQ_SEED_COMP[2], b),
            )
        )
    return verts


def _faq_midpoint(
    verts: list[tuple[float, float, float]],
    cache: dict[tuple[int, int], int],
    a: int,
    b: int,
) -> int:
    """Encounter-order edge-midpoint minting, mirroring the imperative
    ``midpoint`` EXACTLY in cache keying (sorted pair ``a < b``) and numbering
    (next index on first encounter), but with the midpoint COORDINATE computed via
    the midpoint normalize FAQ (``eval_coeff`` on ``_FAQ_MID_COMP``) instead of the
    inline ``(va + vb) / n``. Bit-identical to ``midpoint`` for ``float64``.
    """
    key = (a, b) if a < b else (b, a)
    idx = cache.get(key)
    if idx is not None:
        return idx
    va = verts[a]
    vb = verts[b]
    bnd = {
        "a1": va[0], "a2": va[1], "a3": va[2],
        "b1": vb[0], "b2": vb[1], "b3": vb[2],
    }
    verts.append(
        (
            eval_coeff(_FAQ_MID_COMP[0], bnd),
            eval_coeff(_FAQ_MID_COMP[1], bnd),
            eval_coeff(_FAQ_MID_COMP[2], bnd),
        )
    )
    idx = len(verts) - 1
    cache[key] = idx
    return idx


def duo_subdivide_faq(level: int) -> tuple[np.ndarray, np.ndarray]:
    """Build the level-``level`` DUO icosahedral mesh on the UNIT sphere.

    Returns ``(V (3, Nv) float64, F (3, Nc) int64)`` — a drop-in replacement for
    the imperative ``_subdivide_icosahedron(level)`` (:mod:`duo`). Every
    coordinate is computed through the ESS front-door evaluator ``eval_coeff``
    (the seed normalize FAQ for the 12 base vertices, the midpoint normalize FAQ
    for every refined vertex), but the vertex numbering and face order are the
    imperative ENCOUNTER order — NOT the canonical skolem/rank order. The result
    is byte-identical (``float64``) and same-order as
    ``_subdivide_icosahedron(level)``, column-for-column, value-for-value.
    """
    # Pure-integer structural face table reused from :mod:`duo` (no floats, so the
    # reuse is exact — no transcription risk). Lazy import avoids an import cycle.
    from .duo import _icosahedron_faces

    verts = _faq_seed_vertices()
    F0 = _icosahedron_faces()
    faces: list[tuple[int, int, int]] = [
        (int(F0[0, c]), int(F0[1, c]), int(F0[2, c])) for c in range(F0.shape[1])
    ]

    # Refine: replicate ``_subdivide_icosahedron``'s loop EXACTLY — fresh dict
    # cache per level, iterate faces in order, mint ab/bc/ca via the
    # encounter-order midpoint FAQ, push the 4 children in the imperative order.
    for _ in range(level):
        cache: dict[tuple[int, int], int] = {}
        new_faces: list[tuple[int, int, int]] = []
        for a, b, c in faces:
            ab = _faq_midpoint(verts, cache, a, b)
            bc = _faq_midpoint(verts, cache, b, c)
            ca = _faq_midpoint(verts, cache, c, a)
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
