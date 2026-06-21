"""
    edge_dual_geometry_faq.jl — DUO edge + dual GEOMETRY via the landed
    `EarthSciSerialization` AST evaluator (`eval_coeff`, src/rule_eval.jl).

Declarative companion: `discretizations/grids/duo/faq/edge_dual_geometry.esm` (D2b /
esd-heg.7). That document expresses the DUO Tier-U edge geometry and the spherical
circumcenters of `build_duo_grid` (src/grids/duo.jl) — the cross-product circumcenter
direction (a×b + b×c + c×a, normalized), the per-edge great-circle arcs `R·acos(·)`
of `edge_length` / `cell_distance`, and the geographic `atan2`/`asin` of the
circumcenter — as the RFC `semiring-faq-unified-ir` §8.1 geometry-FAQ pattern.

This file is the thin **evaluation bridge** for the DUO pieces of D2b that
`build_duo_grid` and the Tier-U accessors need. It routes every piece of grid
ARITHMETIC through ESS's single AST pathway (`eval_coeff`); ESD hosts no shadow
evaluator (AGENTS.md single-pathway invariant), so the determinism / cross-binding
byte-identity contract lives entirely in `EarthSciSerialization`. The only ESD-side
logic is the integer connectivity, the circumcenter orientation **sign** decision,
and the boundary-sentinel branch — exactly the parts that are not float arithmetic.

Byte-identity notes (T=Float64):
  * Squares are products (`x*x`, not `^`) — bit-identical to Julia's `literal_pow`.
  * Dot products / the cross-product sum fold in space order, matching the imperative.
  * The arc inputs are pre-divided by R (`P[k,i]/R`) before the dot, then the AST
    applies `R·acos(dot)` — exactly `_duo_arc` (src/grids/duo.jl).
  * The circumcenter's universal outward-orientation flip (`cc·centroid < 0 → negate`)
    is reproduced *conditionally* in plain Julia — the dot is computed in Julia and
    unary `-` is applied to the eval_coeff'd cross-sum components — matching
    `build_duo_grid`'s exact operation order. Unary `-x` is bit-identical to the
    imperative negation, INCLUDING the `−0.0` a negated zero component carries, so the
    normalized direction and its `atan2`/`asin` image are byte-identical.
  * The `clamp` / `max` guards are dropped: on a valid icosahedral mesh every `acos`
    argument already lies in range, so the guards never fire and the bytes match
    (proven bit-for-bit in test/test_duo_edge_dual_geometry_faq.jl, levels 1–3).
"""

_edg_mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

# --- Standard (non-negated) circumcenter cross sum a×b + b×c + c×a, folded in the
#     exact (t1)+(t2)+(t3) order of build_duo_grid. Corner leaves: a=corner1,
#     b=corner2, c=corner3; component k: 1=x,2=y,3=z. The orientation flip is applied
#     CONDITIONALLY in Julia (not baked into the AST) to match the imperative path.
const _EDG_CCX = _edg_mk("+",
    _edg_mk("-", _edg_mk("*", "a2", "b3"), _edg_mk("*", "a3", "b2")),
    _edg_mk("-", _edg_mk("*", "b2", "c3"), _edg_mk("*", "b3", "c2")),
    _edg_mk("-", _edg_mk("*", "c2", "a3"), _edg_mk("*", "c3", "a2")))
const _EDG_CCY = _edg_mk("+",
    _edg_mk("-", _edg_mk("*", "a3", "b1"), _edg_mk("*", "a1", "b3")),
    _edg_mk("-", _edg_mk("*", "b3", "c1"), _edg_mk("*", "b1", "c3")),
    _edg_mk("-", _edg_mk("*", "c3", "a1"), _edg_mk("*", "c1", "a3")))
const _EDG_CCZ = _edg_mk("+",
    _edg_mk("-", _edg_mk("*", "a1", "b2"), _edg_mk("*", "a2", "b1")),
    _edg_mk("-", _edg_mk("*", "b1", "c2"), _edg_mk("*", "b2", "c1")),
    _edg_mk("-", _edg_mk("*", "c1", "a2"), _edg_mk("*", "c2", "a1")))

# --- Norm of the (post-flip) circumcenter direction and the geographic image of the
#     normalized direction. Leaves: the cross-sum components ccx/ccy/ccz and n.
#     Mirrors `n = sqrt(ccx^2+ccy^2+ccz^2)` then `_cart_to_lonlat(ccx/n,ccy/n,ccz/n)`.
const _EDG_CCN = _edg_mk("sqrt",
    _edg_mk("+", _edg_mk("*", "ccx", "ccx"), _edg_mk("*", "ccy", "ccy"), _edg_mk("*", "ccz", "ccz")))   # squares as *, not ^
const _EDG_CC_LON = _edg_mk("atan2", _edg_mk("/", "ccy", "n"), _edg_mk("/", "ccx", "n"))
const _EDG_CC_LAT = _edg_mk("asin", _edg_mk("/", "ccz", "n"))             # no clamp

# --- great-circle arc R·acos(p·q) over point leaves p1..p3, q1..q3 (no clamp); the
#     points are pre-divided by R, mirroring `_duo_arc` (src/grids/duo.jl).
const _EDG_ARC = _edg_mk("*", "R", _edg_mk("acos",
    _edg_mk("+", _edg_mk("*", "p1", "q1"), _edg_mk("*", "p2", "q2"), _edg_mk("*", "p3", "q3"))))

"""
    duo_circumcenter_geo_faq(::Type{T}, vertices_unit, faces, cell_cart) where {T}
        -> (cc_lon, cc_lat)

Spherical-circumcenter geographic coordinates for each DUO triangular cell, from the
declarative FAQ (`discretizations/grids/duo/faq/edge_dual_geometry.esm`), routing all
float arithmetic through the landed ESS AST evaluator (`eval_coeff`).

Arguments:
- `vertices_unit::AbstractMatrix{<:Real}` — `(3, Nv)` **unit-sphere** vertex coords.
- `faces::AbstractMatrix{<:Integer}` — `(3, Nc)` 1-based vertex indices per cell.
- `cell_cart::AbstractMatrix` — `(3, Nc)` R-scaled cell-center cartesian coords, used
  ONLY for the outward-orientation flip-sign dot (not for any returned value).

Returns a `NamedTuple`:
- `cc_lon :: Vector{T}` `(Nc,)` — circumcenter longitude (rad).
- `cc_lat :: Vector{T}` `(Nc,)` — circumcenter latitude (rad).

Byte-identical (T=Float64) to `build_duo_grid`'s circumcenter loop (src/grids/duo.jl
~378-394). The cross-product sum is evaluated via `eval_coeff`; the flip-sign dot and
the conditional unary negation are done in Julia exactly as the imperative code (so
the `±0.0` of a negated zero component is preserved); the norm and the `atan2`/`asin`
image are evaluated via `eval_coeff`.
"""
function duo_circumcenter_geo_faq(::Type{T}, vertices_unit::AbstractMatrix{<:Real},
        faces::AbstractMatrix{<:Integer}, cell_cart::AbstractMatrix) where {T}
    V = vertices_unit
    F = faces
    Nc = size(F, 2)

    cc_lon = Vector{T}(undef, Nc)
    cc_lat = Vector{T}(undef, Nc)

    @inbounds for c in 1:Nc
        a_i = F[1, c]; b_i = F[2, c]; c_i = F[3, c]
        b = Dict{String, Float64}(
            "a1" => Float64(V[1, a_i]), "a2" => Float64(V[2, a_i]), "a3" => Float64(V[3, a_i]),
            "b1" => Float64(V[1, b_i]), "b2" => Float64(V[2, b_i]), "b3" => Float64(V[3, b_i]),
            "c1" => Float64(V[1, c_i]), "c2" => Float64(V[2, c_i]), "c3" => Float64(V[3, c_i]),
        )
        ccx = eval_coeff(_EDG_CCX, b)
        ccy = eval_coeff(_EDG_CCY, b)
        ccz = eval_coeff(_EDG_CCZ, b)
        # Centroid sanity check — flip if facing wrong way. Sign logic + unary negation
        # in Julia, matching build_duo_grid exactly (preserves −0.0 on a negated zero).
        if ccx * Float64(cell_cart[1, c]) + ccy * Float64(cell_cart[2, c]) + ccz * Float64(cell_cart[3, c]) < 0
            ccx = -ccx; ccy = -ccy; ccz = -ccz
        end
        nb = Dict{String, Float64}("ccx" => ccx, "ccy" => ccy, "ccz" => ccz)
        n = eval_coeff(_EDG_CCN, nb)
        nb["n"] = n
        cc_lon[c] = T(eval_coeff(_EDG_CC_LON, nb))
        cc_lat[c] = T(eval_coeff(_EDG_CC_LAT, nb))
    end

    return (cc_lon = cc_lon, cc_lat = cc_lat)
end

# Great-circle arc (radius R) between columns `i` and `j` of an R-scaled cartesian
# array `P`, via the FAQ arc AST. Points are pre-divided by R (matching `_duo_arc`),
# then the AST applies `R·acos(dot)`. clamp is a noop on a valid mesh (dropped).
@inline function _edg_arc_faq(::Type{T}, P::AbstractMatrix, i::Int, j::Int, Rf::Float64) where {T}
    b = Dict{String, Float64}(
        "p1" => Float64(P[1, i]) / Rf, "p2" => Float64(P[2, i]) / Rf, "p3" => Float64(P[3, i]) / Rf,
        "q1" => Float64(P[1, j]) / Rf, "q2" => Float64(P[2, j]) / Rf, "q3" => Float64(P[3, j]) / Rf,
        "R" => Rf,
    )
    return T(eval_coeff(_EDG_ARC, b))
end

"""
    duo_edge_length_faq(::Type{T}, vertices_scaled, edges, R) where {T} -> Vector{T}

Great-circle arc length of each DUO mesh edge (vertex-to-vertex), `(Ne,)`,
edge-aligned with `edges`, from the declarative FAQ
(`discretizations/grids/duo/faq/edge_dual_geometry.esm`), routing the arc arithmetic
through `eval_coeff`.

Arguments:
- `vertices_scaled::AbstractMatrix` — `(3, Nv)` **R-scaled** vertex coordinates.
- `edges::AbstractMatrix{<:Integer}` — `(2, Ne)` vertex-index pairs per edge.
- `R::Real` — sphere radius.

Byte-identical (T=Float64) to `edge_length(g)` (src/grids/duo.jl): replicates
`_duo_arc` exactly — divide each endpoint by R, dot, then `R·acos(·)`; the `clamp`
is a noop on a valid mesh and is dropped.
"""
function duo_edge_length_faq(::Type{T}, vertices_scaled::AbstractMatrix,
        edges::AbstractMatrix{<:Integer}, R::Real) where {T}
    Ne = size(edges, 2)
    Rf = Float64(R)
    out = Vector{T}(undef, Ne)
    @inbounds for e in 1:Ne
        out[e] = _edg_arc_faq(T, vertices_scaled, edges[1, e], edges[2, e], Rf)
    end
    return out
end

"""
    duo_cell_distance_faq(::Type{T}, cell_cart_scaled, edge_cells, R) where {T} -> Vector{T}

Great-circle distance between the centers of the two cells sharing each DUO edge,
`(Ne,)`, edge-aligned with `edge_cells`, from the declarative FAQ
(`discretizations/grids/duo/faq/edge_dual_geometry.esm`), routing the arc arithmetic
through `eval_coeff`.

Arguments:
- `cell_cart_scaled::AbstractMatrix` — `(3, Nc)` **R-scaled** cell-center coords.
- `edge_cells::AbstractMatrix{<:Integer}` — `(2, Ne)` the two incident cell indices
  per edge (`0` for a boundary edge's missing cell).
- `R::Real` — sphere radius.

Byte-identical (T=Float64) to `cell_distance(g)` (src/grids/duo.jl): per edge `e`,
`c1 = edge_cells[1,e]`, `c2 = edge_cells[2,e]`; if both are non-zero, the `_duo_arc`
between the two cell centers, else `zero(T)`. The sentinel branch is plain Julia; the
arc arithmetic is the FAQ AST.
"""
function duo_cell_distance_faq(::Type{T}, cell_cart_scaled::AbstractMatrix,
        edge_cells::AbstractMatrix{<:Integer}, R::Real) where {T}
    Ne = size(edge_cells, 2)
    Rf = Float64(R)
    out = Vector{T}(undef, Ne)
    @inbounds for e in 1:Ne
        c1 = edge_cells[1, e]; c2 = edge_cells[2, e]
        out[e] = (c1 != 0 && c2 != 0) ? _edg_arc_faq(T, cell_cart_scaled, c1, c2, Rf) : zero(T)
    end
    return out
end
