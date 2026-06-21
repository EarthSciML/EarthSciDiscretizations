"""
MPAS-SCVT spherical-topology LEAF (epic `esd-e5m`, bead `esd-e5m.3` / D3).

The one-time, post-convergence topology leaf of the declarative MPAS-SCVT mesh
generator. Given the converged SCVT generators (the cell centres), it emits the
**spherical Delaunay triangulation** of the generators together with its dual
**Voronoi connectivity** (`cells_on_cell` / `edges_on_cell` / the
circumcentre-ring dual `vertices_on_cell`) — the integer topology of an
MPAS-style Voronoi mesh.

# Why this is an irreducible leaf (not a FAQ)

Per the epic's DECLARATIVE-OR-FAIL contract the per-iteration SCVT *step* (D2) is
a declarative semiring-FAQ, but the topology is "an allowed irreducible leaf
(like `intersect_polygon`)". The boundary is drawn exactly here:

  generators ──[LEAF: spherical Delaunay = convex hull of unit vectors]──▶
      triangles + circumcentres ──[DECLARATIVE FAQ: `voronoi_dual_topology_faq`,
      esd-heg.2]──▶ cells_on_cell / edges_on_cell / dual ring

The triangulation step is the part that genuinely cannot be expressed as a
semiring aggregate: it is an incremental, control-flow-heavy convex-hull
construction whose robustness lives in *how* the orientation arithmetic is
evaluated — the same status as `acos`/`sqrt` or the ESS `intersect_polygon`
geometry leaf. Everything *below* the leaf (the angular ordering + bridge-vertex
set-intersection that turn a triangulation into MPAS connectivity) is already the
landed declarative `voronoi_dual_topology_faq`; this leaf simply supplies the
arbitrary-generator triangulation that the DUO primal previously hard-wired.

# Canonical executor: the S2B FFI

In production the triangulation is emitted by the s2bindings.rs S2B FFI
(`SphericalDelaunay::from_lon_lat` → `voronoi_connectivity`, bead `s2b-s7b`),
which resolves orientation with an exact predicate (Shewchuk static filter →
`ExactFloat`) so the integer connectivity is byte-identical across bindings even
at the cospherical degeneracies (cocircular quads) the contract calls out. This
Julia reference binding computes the same hull with a Float64 orientation
predicate — exact for the non-degenerate inputs the determinism contract requires
(see `discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md`) — and emits the
**identical canonical ordering** (CCW-from-outside, smallest-index-first,
lexicographically sorted), so the two agree byte-for-byte on a well-posed mesh.
The contract is the cross-binding glue, exactly as `intersect_polygon` uses
GeometryOps in Julia and `s2bindings` in Rust under one tolerance contract.

The function is a value-free structural operation: connectivity is a pure
function of the input generator coordinates, with no per-step state — it runs
ONCE at convergence, never in the Lloyd recurrence.
"""

# ---------------------------------------------------------------------------
# 3-D convex hull of points in strictly convex position (all on a sphere about
# the origin), i.e. the spherical Delaunay triangulation of the directions.
#
# Incremental insertion: seed a tetrahedron from the first four non-coplanar
# generators, then add each remaining generator by deleting the faces it can
# "see" (its outward side) and re-triangulating the resulting horizon to it.
# Every generator lies on the sphere, so every generator is a hull vertex and
# the hull is a closed simplicial surface (each edge in exactly two triangles).
#
# Determinism: for generators in general position (no four cospherical-coplanar,
# guaranteed by the contract for a well-posed SCVT mesh) the convex hull is
# UNIQUE, so the emitted triangle SET is independent of insertion order; the
# canonical re-ordering applied afterwards (CCW-from-outside, smallest-index
# first, lexicographic sort) makes the emitted ARRAY independent of order too.
# ---------------------------------------------------------------------------

@inline _v3(U, i) = (Float64(U[1, i]), Float64(U[2, i]), Float64(U[3, i]))
@inline _sub(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline function _cross3(a, b)
    return (
        a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1],
    )
end

# Raw (unoriented) face normal of triangle (i,j,k): (Uj-Ui)×(Uk-Ui).
@inline function _face_normal(U, i, j, k)
    a = _v3(U, i)
    return _cross3(_sub(_v3(U, j), a), _sub(_v3(U, k), a))
end

# Signed height of point p above face (i,j,k)'s plane, with the normal oriented
# OUTWARD relative to a strictly-interior reference point `o`. The face is
# "visible" from p (p is on its outward side) when this is > 0. A fixed interior
# `o` (the seed-tetra centroid) is robust where the sphere centre would not be:
# a transient seed face can pass through the origin (e.g. the octahedron's z=0
# diagonal face), making an origin-based orientation test ambiguous.
@inline function _visible_height(U, i, j, k, p, o)
    a = _v3(U, i)
    n = _cross3(_sub(_v3(U, j), a), _sub(_v3(U, k), a))
    if _dot3(n, _sub(a, o)) < 0          # orient n away from the interior point
        n = (-n[1], -n[2], -n[3])
    end
    return _dot3(n, _sub(_v3(U, p), a))
end

# Orient a triangle CCW as seen from OUTSIDE the sphere for the canonical output:
# a real hull face of a polytope enclosing the origin has its outward normal
# pointing away from the centre, i.e. n·a ≥ 0. (Applied only to final hull faces,
# which all genuinely enclose the origin — unlike transient seed faces.)
@inline function _orient_ccw_outside(U, t::NTuple{3, Int})
    n = _face_normal(U, t[1], t[2], t[3])
    a = _v3(U, t[1])
    return _dot3(n, a) >= 0 ? t : (t[1], t[3], t[2])
end

# Ascending-sorted triple — the order-independent identity of a face during hull
# construction (orientation is imposed only at output).
@inline function _sorted3(t::NTuple{3, Int})
    a, b, c = t
    a > b && ((a, b) = (b, a))
    b > c && ((b, c) = (c, b))
    a > b && ((a, b) = (b, a))
    return (a, b, c)
end

# Rotate a CCW triple so its smallest vertex index is first (orientation-
# preserving). The canonical per-triangle form is "CCW-from-outside, min first".
@inline function _rotate_min_first(t::NTuple{3, Int})
    if t[1] <= t[2] && t[1] <= t[3]
        return t
    elseif t[2] <= t[1] && t[2] <= t[3]
        return (t[2], t[3], t[1])
    else
        return (t[3], t[1], t[2])
    end
end

function _spherical_delaunay_triangles(U::AbstractMatrix{<:Real})
    N = size(U, 2)
    N >= 4 || throw(ArgumentError(
        "scvt topology leaf: spherical Delaunay needs ≥ 4 generators; got $N"))

    # --- seed tetrahedron: first three points + first non-coplanar fourth -----
    n012 = _face_normal(U, 1, 2, 3)
    (abs(n012[1]) + abs(n012[2]) + abs(n012[3])) > 0 || throw(ArgumentError(
        "scvt topology leaf: first three generators are collinear"))
    k4 = 0
    @inbounds for k in 4:N
        if abs(_dot3(n012, _sub(_v3(U, k), _v3(U, 1)))) > 1e-14
            k4 = k
            break
        end
    end
    k4 != 0 || throw(ArgumentError(
        "scvt topology leaf: all generators are coplanar (degenerate, not a closed mesh)"))

    # Strictly-interior reference: the seed-tetra centroid. It lies inside the
    # seed tetra and stays inside as the hull only grows, so it orients every
    # transient face unambiguously throughout construction.
    o = (
        (Float64(U[1, 1]) + U[1, 2] + U[1, 3] + U[1, k4]) / 4,
        (Float64(U[2, 1]) + U[2, 2] + U[2, 3] + U[2, k4]) / 4,
        (Float64(U[3, 1]) + U[3, 2] + U[3, 3] + U[3, k4]) / 4,
    )

    # Faces are stored as ascending-sorted triples (order-independent identity);
    # orientation is imposed only on output.
    faces = Set{NTuple{3, Int}}()
    for t in ((1, 2, 3), (1, 2, k4), (1, 3, k4), (2, 3, k4))
        push!(faces, _sorted3(t))
    end

    # --- incremental insertion of the remaining generators -------------------
    @inbounds for p in 4:N
        p == k4 && continue
        # Faces visible from p (p strictly on their outward side, w.r.t. o).
        visible = NTuple{3, Int}[]
        for f in faces
            if _visible_height(U, f[1], f[2], f[3], p, o) > 1e-12
                push!(visible, f)
            end
        end
        isempty(visible) && continue   # p on/inside the current hull (duplicate)

        # Horizon = unordered edges incident to exactly one visible face. Faces
        # are sorted, so each edge is already in (min,max) order.
        edge_count = Dict{Tuple{Int, Int}, Int}()
        for f in visible
            for key in ((f[1], f[2]), (f[2], f[3]), (f[1], f[3]))
                edge_count[key] = get(edge_count, key, 0) + 1
            end
            delete!(faces, f)
        end
        for (key, cnt) in edge_count
            cnt == 1 || continue       # interior edge of the visible cap
            push!(faces, _sorted3((key[1], key[2], p)))
        end
    end

    # --- canonicalize: CCW-from-outside, min-index first, lexicographic sort --
    out = Vector{NTuple{3, Int}}(undef, length(faces))
    i = 0
    for f in faces
        i += 1
        out[i] = _rotate_min_first(_orient_ccw_outside(U, f))
    end
    sort!(out)
    return out
end

"""
    scvt_spherical_delaunay(generators; R=6.371e6) -> NamedTuple

The MPAS-SCVT spherical-topology **leaf**: the spherical Delaunay triangulation
of the generators and its per-triangle Voronoi vertices (circumcentres).

`generators` is `(3, n_cells)` Cartesian generator coordinates on the sphere of
radius `R` (the converged SCVT cell centres). Returns a `NamedTuple`:

- `faces::Matrix{Int}` `(3, n_tri)` — the Delaunay triangles as 1-based generator
  triples, each CCW as seen from outside the sphere, rotated smallest-index-first
  and the whole list sorted lexicographically (the determinism contract). For
  `n` generators forming a closed mesh `n_tri == 2n - 4` (Euler).
- `circumcenters::Matrix{Float64}` `(3, n_tri)` — the dual Voronoi vertices, the
  per-triangle circumcentre directions scaled to radius `R`. Circumcentre `t` is
  dual to triangle `t`.
- `n_triangles::Int` — `n_tri`.

The integer `faces` are governed by the byte-identical determinism contract; the
floating-point `circumcenters` by the geometry tolerance contract (see
`TOPOLOGY_LEAF_CONTRACT.md`).
"""
function scvt_spherical_delaunay(generators::AbstractMatrix{<:Real}; R::Real = 6.371e6)
    size(generators, 1) == 3 || throw(ArgumentError(
        "scvt_spherical_delaunay: generators must be (3, n_cells); got $(size(generators))"))
    (R > 0 && isfinite(R)) ||
        throw(DomainError(R, "scvt_spherical_delaunay: R must be a positive finite number"))
    N = size(generators, 2)

    # Unit directions (the hull is scale-free; normalize so the orientation
    # arithmetic and circumcentres match the contract's unit-sphere convention).
    U = Matrix{Float64}(undef, 3, N)
    @inbounds for i in 1:N
        x = Float64(generators[1, i])
        y = Float64(generators[2, i])
        z = Float64(generators[3, i])
        nrm = sqrt(x * x + y * y + z * z)
        nrm > 0 || throw(ArgumentError(
            "scvt_spherical_delaunay: generator $i is the zero vector"))
        U[1, i] = x / nrm
        U[2, i] = y / nrm
        U[3, i] = z / nrm
    end

    tris = _spherical_delaunay_triangles(U)
    nt = length(tris)
    nt == 2 * N - 4 || throw(ErrorException(
        "scvt_spherical_delaunay: got $nt triangles, expected 2n-4 = $(2 * N - 4) for a " *
            "closed mesh of $N generators (Euler check) — input is not a valid closed " *
            "spherical mesh (degenerate or non-enclosing generators)"))

    faces = Matrix{Int}(undef, 3, nt)
    circ = Matrix{Float64}(undef, 3, nt)
    R_f = Float64(R)
    @inbounds for (t, tri) in enumerate(tris)
        faces[1, t] = tri[1]
        faces[2, t] = tri[2]
        faces[3, t] = tri[3]
        # Circumcentre direction of a spherical triangle with unit-vector corners
        # a,b,c is (a×b + b×c + c×a); normalise and place on the triangle's side.
        a = _v3(U, tri[1])
        b = _v3(U, tri[2])
        c = _v3(U, tri[3])
        cc = _cross3(a, b)
        cc = (cc[1] + (b[2] * c[3] - b[3] * c[2]), cc[2] + (b[3] * c[1] - b[1] * c[3]),
            cc[3] + (b[1] * c[2] - b[2] * c[1]))
        cc = (cc[1] + (c[2] * a[3] - c[3] * a[2]), cc[2] + (c[3] * a[1] - c[1] * a[3]),
            cc[3] + (c[1] * a[2] - c[2] * a[1]))
        # Same hemisphere as the triangle (outward from the sphere centre).
        if _dot3(cc, (a[1] + b[1] + c[1], a[2] + b[2] + c[2], a[3] + b[3] + c[3])) < 0
            cc = (-cc[1], -cc[2], -cc[3])
        end
        cn = sqrt(cc[1]^2 + cc[2]^2 + cc[3]^2)
        circ[1, t] = R_f * cc[1] / cn
        circ[2, t] = R_f * cc[2] / cn
        circ[3, t] = R_f * cc[3] / cn
    end

    return (faces = faces, circumcenters = circ, n_triangles = nt)
end

"""
    scvt_voronoi_connectivity(generators; R=6.371e6) -> NamedTuple

The full MPAS-SCVT spherical-topology leaf output: the spherical Delaunay
triangulation of the `generators` (the irreducible hull leaf) composed with the
landed declarative dual-topology FAQ (`voronoi_dual_topology_faq`, esd-heg.2) to
emit the **deterministic Voronoi connectivity** of the dual mesh.

`generators` is `(3, n_cells)` Cartesian cell centres on the sphere of radius
`R`. Returns a `NamedTuple` (MPAS field conventions, 1-based with `0` pad):

- `faces`, `circumcenters`, `n_triangles` — the leaf triangulation (see
  [`scvt_spherical_delaunay`](@ref)). The circumcentres are the dual Voronoi
  vertices; `n_triangles == n_vertices` of the dual mesh.
- `n_edges_on_cell::Vector{Int}` `(n_cells,)` — valence of each Voronoi cell.
- `cells_on_cell::Matrix{Int}` `(max_edges, n_cells)` — neighbour generators in
  CCW ring order, `0`-padded past the valence.
- `edges_on_cell::Matrix{Int}` `(max_edges, n_cells)` — the canonical shared-edge
  id in the same ring slot.
- `vertices_on_cell::Matrix{Int}` `(max_edges, n_cells)` — the surrounding dual
  Voronoi vertices (circumcentre / triangle indices) in CCW ring order, index-
  aligned with `cells_on_cell` per the S2B convention: `vertices_on_cell[i, c]`
  is the circumcentre of the Delaunay triangle between neighbours
  `cells_on_cell[i, c]` and `cells_on_cell[(i % k) + 1, c]`.
- `max_edges::Int`.

The integer connectivity is deterministic and cross-binding byte-identical by
the determinism contract; the circumcentre geometry rides the tolerance
contract. Runs ONCE at convergence (never in the Lloyd recurrence).
"""
function scvt_voronoi_connectivity(generators::AbstractMatrix{<:Real}; R::Real = 6.371e6)
    leaf = scvt_spherical_delaunay(generators; R = R)
    Ncells = size(generators, 2)

    # Generators scaled to radius R, matching the `cell_cart`/`vertices` (3,Nv)
    # convention the dual-topology FAQ shares with the imperative builder.
    R_f = Float64(R)
    verts = Matrix{Float64}(undef, 3, Ncells)
    @inbounds for i in 1:Ncells
        x = Float64(generators[1, i]); y = Float64(generators[2, i]); z = Float64(generators[3, i])
        nrm = sqrt(x * x + y * y + z * z)
        verts[1, i] = R_f * x / nrm
        verts[2, i] = R_f * y / nrm
        verts[3, i] = R_f * z / nrm
    end

    # Declarative dual topology (angular ordering + bridge-vertex intersection).
    dual = voronoi_dual_topology_faq(verts, leaf.faces, leaf.circumcenters; R = R_f)

    # vertices_on_cell: the dual Voronoi vertices around each cell are the
    # angularly-ordered incident triangles (circumcentres). The S2B convention
    # (s2b-s7b) index-aligns slot i with the Voronoi vertex BETWEEN neighbours
    # `cells_on_cell[i]` and `cells_on_cell[i+1]`. `cells_on_cell[i]` is the bridge
    # of the consecutive triangle pair (f_i, f_{i+1}) in `sorted_vertex_faces`, so
    # the vertex shared by the edges to neighbours i and i+1 is f_{i+1} — i.e.
    # `sorted_vertex_faces[(i % k) + 1]`. Matching this here makes the Julia
    # reference and the D4 S2B-wrap emit byte-identical `vertices_on_cell`.
    vertices_on_cell = zeros(Int, dual.max_edges, Ncells)
    @inbounds for c in 1:Ncells
        ring = dual.sorted_vertex_faces[c]
        kv = length(ring)
        for i in 1:kv
            vertices_on_cell[i, c] = ring[(i % kv) + 1]
        end
    end

    return (
        faces = leaf.faces,
        circumcenters = leaf.circumcenters,
        n_triangles = leaf.n_triangles,
        n_edges_on_cell = dual.n_edges_on_cell,
        cells_on_cell = dual.cells_on_cell,
        edges_on_cell = dual.edges_on_cell,
        vertices_on_cell = vertices_on_cell,
        max_edges = dual.max_edges,
    )
end
