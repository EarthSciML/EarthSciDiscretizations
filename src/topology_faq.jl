"""
    topology_faq.jl — DUO primal triangular-mesh topology via the landed M3
    value-invention engine (`EarthSciSerialization.Relational`).

Declarative companion: `discretizations/grids/duo/rules/primal_topology.esm`
expresses this chain as a value-invention FAQ document (RFC `semiring-faq-unified-ir`
§7.3: enumerate edges from the face→vertex relation, skolem on the sorted
vertex-pair key, `distinct` dedup, `rank` dense ids, then the inversion joins for
`edges_on_face` / `cell_neighbors` and the group-by for `vertex_faces`).

This module is the thin **evaluation bridge**: it shapes a triangular `faces`
matrix into the relational inputs and routes every dedup / rank / join through
ESS's `Relational` primitives. ESD hosts no shadow relational engine — the
determinism contract (and therefore the cross-binding byte-identity of the output)
lives entirely in `EarthSciSerialization` (`CONFORMANCE_SPEC.md` §5.5; AGENTS.md
single-pathway invariant; GRIDS_API §4.3). The only ESD-side logic is the
*authoring* of which vertex pairs are a triangle's edges — the declarative content
of the FAQ rule itself.
"""

import EarthSciSerialization
const _ESS_RELATIONAL = EarthSciSerialization.Relational

"""
    primal_topology_faq(faces, n_vertices) -> NamedTuple

Materialize DUO primal unstructured connectivity from a `(3, Nc)` triangular
`faces` matrix (1-based vertex ids) by evaluating the value-invention FAQ through
the landed ESS M3 engine. Returns a NamedTuple with the four primal-topology
arrays, matching the `DuoGrid` field conventions:

- `edges::Matrix{Int}` `(2, Ne)` — sorted vertex pairs `(a < b)`, emitted in the
  ESS §5.7 canonical (sorted) order. This is the cross-binding byte-stable edge
  numbering; the imperative `_build_connectivity` emitted the *same set* in
  `Dict`-insertion order (non-canonical, binding-divergent), which value-invention
  replaces with the deterministic sorted order.
- `edges_on_face::Matrix{Int}` `(3, Nc)` — 1-based dense id of edge `k` of face
  `c`, where slot `k` is the edge OPPOSITE vertex `k` (`k=1→(v2,v3)`,
  `k=2→(v3,v1)`, `k=3→(v1,v2)`), the same slot convention as `cell_neighbors`.
- `cell_neighbors::Matrix{Int}` `(3, Nc)` — cell sharing edge `k` of face `c`;
  `0` for a boundary edge (none on a closed icosahedral mesh).
- `vertex_faces::Vector{Vector{Int}}` — sorted face ids incident on each vertex.

The dedup (`distinct`), dense renumbering (`rank`), and connectivity inversion
(`equijoin`) are computed by `EarthSciSerialization.Relational`, so the result is
byte-identical to the Rust / Python bindings' M3 engine by ESS's cross-binding
determinism contract.
"""
function primal_topology_faq(faces::AbstractMatrix{<:Integer}, n_vertices::Integer)
    R = _ESS_RELATIONAL
    size(faces, 1) == 3 ||
        throw(ArgumentError("primal_topology_faq: faces must be (3, Nc) triangular; got $(size(faces))"))
    n_vertices >= 0 ||
        throw(DomainError(n_vertices, "primal_topology_faq: n_vertices must be ≥ 0"))
    Nc = size(faces, 2)

    # The three undirected edges of triangle `c`, indexed by the local slot `k`
    # OPPOSITE vertex `k` (k=1→(v2,v3), k=2→(v3,v1), k=3→(v1,v2)) — the slot
    # convention shared by `cell_neighbors`/`edges_on_face` and matched to the
    # imperative `_build_connectivity`.
    _face_edges(c) = let v1 = Int(faces[1, c]), v2 = Int(faces[2, c]), v3 = Int(faces[3, c])
        ((v2, v3), (v3, v1), (v1, v2))
    end

    # --- (1) edges: skolem on the sorted vertex pair + distinct + rank ----------
    # One canonical undirected key per (face, local-edge); `distinct` dedups the
    # whole multiset into the sorted edge set, `rank` assigns dense 1-based ids.
    edge_keys = Tuple{Int, Int}[]
    sizehint!(edge_keys, 3 * Nc)
    for c in 1:Nc, (a, b) in _face_edges(c)
        push!(edge_keys, R.skolem_edge(a, b))
    end
    edge_order = R.distinct(edge_keys)        # sorted unique (min, max) pairs
    ranking = R.rank(edge_keys; base = 1)     # 1-based dense ids (Julia native)
    Ne = length(edge_order)
    edges = Matrix{Int}(undef, 2, Ne)
    @inbounds for (e, ab) in enumerate(edge_order)
        edges[1, e] = ab[1]
        edges[2, e] = ab[2]
    end

    # --- (2) edges_on_face: dense edge id of each of a face's three edges -------
    # `rank` materialized the (key → dense id) map; gathering it per (face, k) is
    # the face→edge value-equality join against that dense numbering.
    edges_on_face = Matrix{Int}(undef, 3, Nc)
    @inbounds for c in 1:Nc
        for (k, ab) in enumerate(_face_edges(c))
            edges_on_face[k, c] = ranking.id[R.skolem_edge(ab[1], ab[2])]
        end
    end

    # --- (3) cell_neighbors: the OTHER cell across each edge --------------------
    # Inversion join (ESS `equijoin` — the documented connectivity-inversion
    # primitive): self-join the (dense edge id, cell, slot) incidence on the edge
    # id; each pair `(c1, c2)` with `c1 ≠ c2` on a shared edge is the neighbour
    # across that slot. A closed icosahedral mesh is manifold (2 cells per edge).
    incidence = Tuple{Int, Int, Int}[]       # (dense edge id, cell, slot k)
    sizehint!(incidence, 3 * Nc)
    @inbounds for c in 1:Nc
        for (k, ab) in enumerate(_face_edges(c))
            push!(incidence, (ranking.id[R.skolem_edge(ab[1], ab[2])], c, k))
        end
    end
    cell_neighbors = zeros(Int, 3, Nc)
    joined = R.equijoin(incidence, incidence; on_left = t -> t[1], on_right = t -> t[1])
    @inbounds for pair in joined
        l, r = pair[1], pair[2]
        c1, k1 = l[2], l[3]
        c2 = r[2]
        c1 == c2 && continue
        if cell_neighbors[k1, c1] != 0 && cell_neighbors[k1, c1] != c2
            throw(AssertionError(
                "primal_topology_faq: non-manifold edge — cell $c1 slot $k1 shares an edge with >2 cells"))
        end
        cell_neighbors[k1, c1] = c2
    end

    # --- (4) vertex_faces: faces incident on each vertex, sorted ----------------
    # Group-by over the face→vertex relation: `distinct` over (vertex, face) rows
    # sorts by (vertex, face), so each vertex's faces land contiguous and ascending
    # — the ragged faces-incident-to-vertex set.
    vf_rows = Tuple{Int, Int}[]
    sizehint!(vf_rows, 3 * Nc)
    @inbounds for c in 1:Nc, k in 1:3
        push!(vf_rows, (Int(faces[k, c]), c))
    end
    vf_sorted = R.distinct(vf_rows)
    vertex_faces = [Int[] for _ in 1:n_vertices]
    @inbounds for vf in vf_sorted
        push!(vertex_faces[vf[1]], vf[2])
    end

    return (
        edges = edges,
        edges_on_face = edges_on_face,
        cell_neighbors = cell_neighbors,
        vertex_faces = vertex_faces,
    )
end
