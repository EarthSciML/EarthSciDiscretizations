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

# Shared empty-bridge sentinel: the `get` fallback for a consecutive face pair
# with no recorded shared vertex (a degenerate fan; the bridge loop then errors).
const EMPTY_BRIDGE = Int[]

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

"""
    voronoi_dual_topology_faq(vertices, faces, cell_cart; R=6.371e6) -> NamedTuple

Materialize the **MPAS Voronoi-dual topology** of a DUO icosahedral primal mesh
as a value-invention FAQ evaluated through the landed ESS M3 engine. Each MPAS
dual cell `v` IS the primal vertex `v`; its neighbours / incident edges are the
ring of primal faces and bridge vertices around `v`. Replaces the topology half
of the imperative `_duo_voronoi_dual` (`src/grids/mpas.jl`, steps 3–4). Returns a
NamedTuple, matching the `MpasMeshData` field conventions:

- `n_edges_on_cell::Vector{Int}` `(Nv,)` — valence of each dual cell = the number
  of primal faces incident on primal vertex `v` (5 at the 12 pentagonal seeds, 6
  elsewhere on a closed icosahedral dual).
- `cells_on_cell::Matrix{Int}` `(max_edges, Nv)` — dual-cell neighbours (the
  primal **bridge vertices**) in angular order, `0`-padded past the valence.
- `edges_on_cell::Matrix{Int}` `(max_edges, Nv)` — the shared primal edge id of
  each neighbour in the SAME angular slot, `0`-padded. Ids are the D1a
  cross-binding **canonical** edge numbering (`primal_topology_faq(...).edges`),
  not the imperative `_build_connectivity` Dict-insertion numbering.
- `max_edges::Int` — `maximum(n_edges_on_cell)` (6 for a closed icosahedral dual).
- `sorted_vertex_faces::Vector{Vector{Int}}` — the angularly-ordered incident
  faces per vertex (the dual-cell corner ring), exposed because the dual geometry
  FAQ (dual-cell areas over the circumcenter ring, bead esd-heg.7 / D2b) consumes
  this exact ordering.

# Two value-invention steps (RFC `semiring-faq-unified-ir` §7.3)

1. **Angular ordering of the incident faces** — a per-vertex *rank* of the primal
   faces incident on `v` (`primal_topology_faq(...).vertex_faces`, the D1a
   group-by) by their tangent-plane `atan2` angle: project each incident face
   centroid (`cell_cart`) into the orthonormal tangent basis at the unit vertex
   direction and sort by `atan2`. This is the one GEOMETRIC step — the angle is a
   numeric `atan2`, NOT a relational key (floats are forbidden in keys, §5.5.1
   rule 1), so the ordering is a stable sort here in the bridge rather than
   `Relational.rank`. It reproduces the imperative `_duo_voronoi_dual` sort
   byte-for-byte.
2. **Bridge-vertex set-INTERSECTION** — for each consecutive face pair
   `(f_i, f_{i+1})` around `v` the dual neighbour is the *bridge vertex* `w`: the
   unique vertex shared by both faces other than `v`. This is the relational
   intersection `V(f_i) ∩ V(f_{i+1}) ∖ {v}`, materialised by `Relational.equijoin`
   of the (vertex, face) incidence on the shared vertex id. The incident edge is
   then `edges_on_cell = rank(skolem_edge(v, w))` — the D1a canonical edge id of
   the primal edge `{v, w}`.

The intersection and the edge-id lookup key on integer ids, so they are
byte-identical across the Julia / Rust / Python bindings by ESS's §5.5
determinism contract; the angular `atan2` sort is pinned per icosahedral level in
`tests/conformance/grids/duo/voronoi_dual_topology/`. Declarative companion:
`discretizations/grids/duo/rules/voronoi_dual_topology.esm`.
"""
function voronoi_dual_topology_faq(
        vertices::AbstractMatrix{<:Real},
        faces::AbstractMatrix{<:Integer},
        cell_cart::AbstractMatrix{<:Real};
        R::Real = 6.371e6,
    )
    Rel = _ESS_RELATIONAL
    size(vertices, 1) == 3 || throw(ArgumentError(
        "voronoi_dual_topology_faq: vertices must be (3, Nv); got $(size(vertices))"))
    size(faces, 1) == 3 || throw(ArgumentError(
        "voronoi_dual_topology_faq: faces must be (3, Nc) triangular; got $(size(faces))"))
    size(cell_cart, 1) == 3 || throw(ArgumentError(
        "voronoi_dual_topology_faq: cell_cart must be (3, Nc); got $(size(cell_cart))"))
    size(cell_cart, 2) == size(faces, 2) || throw(ArgumentError(
        "voronoi_dual_topology_faq: cell_cart $(size(cell_cart)) and faces $(size(faces)) must share Nc"))
    (R > 0 && isfinite(R)) ||
        throw(DomainError(R, "voronoi_dual_topology_faq: R must be a positive finite number"))

    Nv = size(vertices, 2)
    Nc = size(faces, 2)
    R_f = Float64(R)

    # --- D1a reuse: vertex→face group-by + the canonical edge numbering --------
    # The dual topology rides D1a's primal value-invention: `vertex_faces` is the
    # group-by basis the angular sort orders; `edges` is the canonical (sorted)
    # edge set whose dense ids `edges_on_cell` references.
    prim = primal_topology_faq(faces, Nv)
    vertex_faces = prim.vertex_faces
    Ne = size(prim.edges, 2)
    edge_id = Dict{Tuple{Int, Int}, Int}()      # canonical (min,max) primal edge → dense id
    sizehint!(edge_id, Ne)
    @inbounds for e in 1:Ne
        edge_id[(prim.edges[1, e], prim.edges[2, e])] = e
    end

    # --- (1) angular ordering of incident faces (the one geometric step) -------
    sorted_vf = _angular_order_faces(vertices, cell_cart, vertex_faces, R_f)
    n_edges_on_cell = Int[length(sorted_vf[v]) for v in 1:Nv]
    max_edges = isempty(n_edges_on_cell) ? 0 : maximum(n_edges_on_cell)

    # --- (2) bridge-vertex set-intersection via relational equijoin ------------
    # Self-equijoin the (vertex, face) incidence on the shared vertex id: a matched
    # pair (l, r) with faces fa < fb means vertex `l[1]` is shared by fa AND fb.
    # Two faces consecutive around `v` share exactly the primal edge {v, w}, so the
    # shared-vertex set of (f_i, f_{i+1}) is {v, w}; the bridge `w` is the non-`v`
    # member. The join keys on integer vertex ids — ESS guarantees its order is
    # binding-stable, but here only set membership is used.
    fv = Tuple{Int, Int}[]                       # (vertex, face) incidence
    sizehint!(fv, 3 * Nc)
    @inbounds for f in 1:Nc, k in 1:3
        push!(fv, (Int(faces[k, f]), f))
    end
    joined = Rel.equijoin(fv, fv; on_left = t -> t[1], on_right = t -> t[1])
    shared = Dict{Tuple{Int, Int}, Vector{Int}}()
    @inbounds for pair in joined
        l, r = pair[1], pair[2]
        fa, fb = l[2], r[2]
        fa < fb || continue                      # each unordered face pair once per shared vertex
        push!(get!(() -> Int[], shared, (fa, fb)), l[1])
    end

    cells_on_cell = zeros(Int, max_edges, Nv)
    edges_on_cell = zeros(Int, max_edges, Nv)
    @inbounds for v in 1:Nv
        sfaces = sorted_vf[v]
        kv = length(sfaces)
        for i in 1:kv
            fi = sfaces[i]
            fnext = sfaces[(i % kv) + 1]
            key = fi < fnext ? (fi, fnext) : (fnext, fi)
            cand = get(shared, key, EMPTY_BRIDGE)
            w = 0
            nbridge = 0
            for x in cand
                if x != v
                    w = x
                    nbridge += 1
                end
            end
            nbridge == 1 || error(
                "voronoi_dual_topology_faq: vertex $v: expected 1 bridge vertex between " *
                    "consecutive faces $fi/$fnext, got $nbridge")
            cells_on_cell[i, v] = w
            edges_on_cell[i, v] = edge_id[v < w ? (v, w) : (w, v)]
        end
    end

    return (
        n_edges_on_cell = n_edges_on_cell,
        cells_on_cell = cells_on_cell,
        edges_on_cell = edges_on_cell,
        max_edges = max_edges,
        sorted_vertex_faces = sorted_vf,
    )
end

# Angular (CCW-from-outside-the-sphere) ordering of the faces incident on each
# vertex, by the tangent-plane `atan2` angle of the face centroids — the one
# geometric step of the dual-topology value-invention. Transcribes the imperative
# `_duo_voronoi_dual` sort (`src/grids/mpas.jl` step 3) operation-for-operation so
# the resulting `cells_on_cell` order is byte-identical. A stable sort by a numeric
# `atan2` angle is NOT a relational key op (floats are forbidden in keys, §5.5.1
# rule 1); its cross-binding determinism comes from the geometry pipeline's fixed
# evaluation order, pinned by the conformance golden.
function _angular_order_faces(
        vertices::AbstractMatrix{<:Real},
        cell_cart::AbstractMatrix{<:Real},
        vertex_faces::Vector{Vector{Int}},
        R_f::Float64,
    )
    Nv = length(vertex_faces)
    sorted_vf = Vector{Vector{Int}}(undef, Nv)
    @inbounds for v in 1:Nv
        faces_v = vertex_faces[v]
        kv = length(faces_v)
        if kv <= 1
            sorted_vf[v] = copy(faces_v)
            continue
        end
        Pvx = Float64(vertices[1, v]) / R_f
        Pvy = Float64(vertices[2, v]) / R_f
        Pvz = Float64(vertices[3, v]) / R_f
        # Tangent-plane orthonormal basis (e1, e2) at the unit vertex direction.
        if abs(Pvz) < 0.9
            d = Pvz                               # ref = ẑ; (ref·Pv) = Pvz
            e1x, e1y, e1z = -d * Pvx, -d * Pvy, 1.0 - d * Pvz
        else
            d = Pvx                               # ref = x̂
            e1x, e1y, e1z = 1.0 - d * Pvx, -d * Pvy, -d * Pvz
        end
        n1 = sqrt(e1x^2 + e1y^2 + e1z^2)
        e1x /= n1; e1y /= n1; e1z /= n1
        e2x = Pvy * e1z - Pvz * e1y               # e2 = Pv × e1 (unit, Pv ⊥ e1)
        e2y = Pvz * e1x - Pvx * e1z
        e2z = Pvx * e1y - Pvy * e1x

        angles = Vector{Float64}(undef, kv)
        for (i, f) in enumerate(faces_v)
            qx = Float64(cell_cart[1, f]) / R_f
            qy = Float64(cell_cart[2, f]) / R_f
            qz = Float64(cell_cart[3, f]) / R_f
            dqPv = qx * Pvx + qy * Pvy + qz * Pvz
            qpx = qx - dqPv * Pvx                 # project the centroid into the tangent plane
            qpy = qy - dqPv * Pvy
            qpz = qz - dqPv * Pvz
            angles[i] = atan(
                qpx * e2x + qpy * e2y + qpz * e2z,
                qpx * e1x + qpy * e1y + qpz * e1z,
            )
        end
        sorted_vf[v] = faces_v[sortperm(angles)]
    end
    return sorted_vf
end
