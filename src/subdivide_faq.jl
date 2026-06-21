"""
    subdivide_faq.jl — DUO icosahedral mesh construction with every vertex
    COORDINATE routed through the ESS front-door evaluator (`eval_coeff`),
    while preserving the EXACT combinatorial ORDER of the imperative
    `_subdivide_icosahedron` (src/grids/duo.jl).

Declarative companions: `discretizations/grids/duo/faq/subdivide_fold.esm` and
`discretizations/grids/duo/faq/subdivide_refine.esm` (D3). Those documents
express the level-fold and the one-level value-invention refine pass as
semiring-FAQ models; the seed `vert_raw / sqrt(seed_sq)` normalize and the
midpoint `(a+b) / sqrt(Σ (a+b)²)` average are the geometry FAQs this file
evaluates through `eval_coeff` (src/rule_eval.jl — ESD's single arithmetic
pathway; no shadow evaluator, per AGENTS.md / GRIDS_API §4.3).

## Front-door coordinates, imperative encounter-order numbering

This routine routes only the **coordinates** through the front door: every
vertex value (the 12 seed vertices and every refined midpoint) is computed by
`eval_coeff` on the SAME ASTs the landed conformance tests drive
(`test/test_duo_subdivision_fold.jl` `_rcomp`, `test/test_duo_subdivision_faq.jl`
`_component`). The **numbering**, however, is the imperative encounter order —
12 base vertices in table order, then midpoints minted in `Dict`-cache
ENCOUNTER order as faces are iterated — NOT the canonical sorted/skolem/rank
order serialized in `subdivide_refine.esm`.

The canonical skolem/rank order is the **test-only cross-binding artifact**: it
reorders the invented midpoints into a binding-neutral sequence for byte-golden
comparison. The PRODUCTION path must keep the imperative order because
downstream dual-conformance (`test/test_duo_voronoi_dual_topology_faq.jl:58-59`)
asserts the per-vertex dual arrays equal the imperative `_duo_voronoi_dual`
EXACTLY, and the dual is indexed by primal vertex — so the primal vertex
numbering must match `_subdivide_icosahedron` column-for-column. `eval_coeff` on
these ASTs is bit-identical to the imperative `_icosahedron_vertices` /
`_midpoint!` arithmetic — Julia lowers `x^2` to `x*x`, so `(a+b)*(a+b)` equals
the imperative `(a+b)^2` to the ULP, the 3-arg `+` folds left-assoc exactly as
`v1²+v2²+v3²`, and the single divide-by-`sqrt` is the same op — so the output is
byte-identical AND same-order as `_subdivide_icosahedron(T, level)` for T=Float64.
"""

# Raw golden-ratio icosahedron vertices, pre-normalization. Transcribed to match
# `_icosahedron_vertices`'s `raw` literal table (src/grids/duo.jl) row-for-row;
# φ = (1 + sqrt(5))/2. Stored as a module const so it is not rebuilt per call.
const _FAQ_PHI = (1.0 + sqrt(5.0)) / 2.0
const _FAQ_VERT_RAW = [
    0.0  1.0       _FAQ_PHI;
    0.0 -1.0       _FAQ_PHI;
    0.0  1.0      -_FAQ_PHI;
    0.0 -1.0      -_FAQ_PHI;
    1.0   _FAQ_PHI 0.0;
    -1.0  _FAQ_PHI 0.0;
    1.0  -_FAQ_PHI 0.0;
    -1.0 -_FAQ_PHI 0.0;
    _FAQ_PHI 0.0   1.0;
    _FAQ_PHI 0.0  -1.0;
    -_FAQ_PHI 0.0  1.0;
    -_FAQ_PHI 0.0 -1.0;
]

# --- Seed normalize-to-sphere FAQ AST -----------------------------------------
# coord_d = r$d / sqrt(r1*r1 + r2*r2 + r3*r3).  Squares are products (`*`), NOT
# `^`, and the 3-arg `+` folds left-associatively — the exact float-op sequence
# of the imperative `n = sqrt(v[1]^2+v[2]^2+v[3]^2); v[d]/n` (Julia `v^2` lowers
# to `v*v`). Held as module consts (the same shape as `test/test_duo_subdivision_fold.jl`
# `_rcomp`) so the per-vertex evaluation does not re-allocate the trees.
_faq_rsq(n) = Dict{String, Any}("op" => "*", "args" => Any[n, n])
const _FAQ_SEED_NORM = Dict{String, Any}(
    "op" => "sqrt",
    "args" => Any[
        Dict{String, Any}(
            "op" => "+",
            "args" => Any[_faq_rsq("r1"), _faq_rsq("r2"), _faq_rsq("r3")],
        ),
    ],
)
const _FAQ_SEED_COMP = (
    Dict{String, Any}("op" => "/", "args" => Any["r1", _FAQ_SEED_NORM]),
    Dict{String, Any}("op" => "/", "args" => Any["r2", _FAQ_SEED_NORM]),
    Dict{String, Any}("op" => "/", "args" => Any["r3", _FAQ_SEED_NORM]),
)

# --- Midpoint normalize-to-sphere FAQ AST -------------------------------------
# coord_d = (a$d + b$d) / sqrt(Σ_d (a$d + b$d)*(a$d + b$d)).  Again squares are
# products so the float ops match the imperative `_midpoint!` bit-for-bit
# (`mx^2` lowers to `mx*mx`). Same shape as `test/test_duo_subdivision_faq.jl`
# `_component`.
_faq_sum(d) = Dict{String, Any}("op" => "+", "args" => Any["a$d", "b$d"])
_faq_sq(d) = Dict{String, Any}("op" => "*", "args" => Any[_faq_sum(d), _faq_sum(d)])
const _FAQ_MID_NORM = Dict{String, Any}(
    "op" => "sqrt",
    "args" => Any[
        Dict{String, Any}(
            "op" => "+",
            "args" => Any[_faq_sq(1), _faq_sq(2), _faq_sq(3)],
        ),
    ],
)
const _FAQ_MID_COMP = (
    Dict{String, Any}("op" => "/", "args" => Any[_faq_sum(1), _FAQ_MID_NORM]),
    Dict{String, Any}("op" => "/", "args" => Any[_faq_sum(2), _FAQ_MID_NORM]),
    Dict{String, Any}("op" => "/", "args" => Any[_faq_sum(3), _FAQ_MID_NORM]),
)

"""
    _faq_seed_vertices(::Type{T}) -> Vector{NTuple{3,T}}

The 12 base icosahedron vertices on the unit sphere, in `_FAQ_VERT_RAW` table
order (= `_icosahedron_vertices` order), each computed via the seed normalize
FAQ (`eval_coeff`). Bit-identical to `_icosahedron_vertices(T)` columns for
T=Float64.
"""
function _faq_seed_vertices(::Type{T}) where {T}
    Nv = size(_FAQ_VERT_RAW, 1)
    verts = Vector{NTuple{3, T}}(undef, Nv)
    @inbounds for i in 1:Nv
        b = Dict{String, Float64}(
            "r1" => Float64(_FAQ_VERT_RAW[i, 1]),
            "r2" => Float64(_FAQ_VERT_RAW[i, 2]),
            "r3" => Float64(_FAQ_VERT_RAW[i, 3]),
        )
        verts[i] = (
            T(eval_coeff(_FAQ_SEED_COMP[1], b)),
            T(eval_coeff(_FAQ_SEED_COMP[2], b)),
            T(eval_coeff(_FAQ_SEED_COMP[3], b)),
        )
    end
    return verts
end

"""
    _faq_midpoint!(verts, cache, a, b) -> Int

Encounter-order edge-midpoint minting, mirroring `_midpoint!` (src/grids/duo.jl)
EXACTLY in cache keying (sorted pair `a<b`) and numbering (next index on first
encounter), but with the midpoint COORDINATE computed via the midpoint
normalize FAQ (`eval_coeff` on `_FAQ_MID_COMP`) instead of the inline
`(va+vb)/n`. Bit-identical to `_midpoint!` for T=Float64.
"""
function _faq_midpoint!(
        verts::Vector{NTuple{3, T}}, cache::Dict{Tuple{Int, Int}, Int},
        a::Int, b::Int
    ) where {T}
    key = a < b ? (a, b) : (b, a)
    idx = get(cache, key, 0)
    if idx != 0
        return idx
    end
    va = verts[a]
    vb = verts[b]
    bnd = Dict{String, Float64}(
        "a1" => Float64(va[1]), "a2" => Float64(va[2]), "a3" => Float64(va[3]),
        "b1" => Float64(vb[1]), "b2" => Float64(vb[2]), "b3" => Float64(vb[3]),
    )
    push!(
        verts,
        (
            T(eval_coeff(_FAQ_MID_COMP[1], bnd)),
            T(eval_coeff(_FAQ_MID_COMP[2], bnd)),
            T(eval_coeff(_FAQ_MID_COMP[3], bnd)),
        ),
    )
    idx = length(verts)
    cache[key] = idx
    return idx
end

"""
    duo_subdivide_faq(::Type{T}, level::Int) where {T} -> (V::Matrix{T}, F::Matrix{Int})

Build the level-`level` DUO icosahedral mesh on the UNIT sphere, returning
vertices `V` (`3×Nv`, element type `T`) and faces `F` (`3×Nc`, `Int`).

Every coordinate is computed through the ESS front-door evaluator `eval_coeff`
(the seed normalize FAQ for the 12 base vertices, the midpoint normalize FAQ for
every refined vertex), but the vertex numbering and face order are the imperative
ENCOUNTER order of `_subdivide_icosahedron` (src/grids/duo.jl) — NOT the
canonical skolem/rank order. The result is byte-identical (for T=Float64) and
same-order as `_subdivide_icosahedron(T, level)`, column-for-column,
value-for-value. See the file docstring for why downstream dual byte-identity
requires the imperative order.
"""
function duo_subdivide_faq(::Type{T}, level::Int) where {T}
    # Seed: 12 base vertices via the seed FAQ (table order); 20 base faces from
    # the pure-integer structural table (`_icosahedron_faces` carries no floats,
    # so reusing it is exact — no transcription risk).
    verts = _faq_seed_vertices(T)
    F0 = _icosahedron_faces()
    faces = Tuple{Int, Int, Int}[(F0[1, c], F0[2, c], F0[3, c]) for c in 1:size(F0, 2)]

    # Refine: replicate `_subdivide_icosahedron`'s loop EXACTLY — fresh Dict
    # cache per level, iterate faces in order, mint ab/bc/ca via the encounter-
    # order midpoint FAQ, push the 4 children in the imperative order.
    for _ in 1:level
        cache = Dict{Tuple{Int, Int}, Int}()
        new_faces = Tuple{Int, Int, Int}[]
        sizehint!(new_faces, 4 * length(faces))
        for (a, b, c) in faces
            ab = _faq_midpoint!(verts, cache, a, b)
            bc = _faq_midpoint!(verts, cache, b, c)
            ca = _faq_midpoint!(verts, cache, c, a)
            push!(new_faces, (a, ab, ca))
            push!(new_faces, (b, bc, ab))
            push!(new_faces, (c, ca, bc))
            push!(new_faces, (ab, bc, ca))
        end
        faces = new_faces
    end

    Nv = length(verts)
    Nc = length(faces)
    V = Matrix{T}(undef, 3, Nv)
    for i in 1:Nv
        v = verts[i]
        V[1, i] = v[1]
        V[2, i] = v[2]
        V[3, i] = v[3]
    end
    F = Matrix{Int}(undef, 3, Nc)
    for c in 1:Nc
        f = faces[c]
        F[1, c] = f[1]
        F[2, c] = f[2]
        F[3, c] = f[3]
    end
    return V, F
end
