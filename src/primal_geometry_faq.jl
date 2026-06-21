"""
    primal_geometry_faq.jl — DUO icosahedral primal-cell geometry via the landed
    `EarthSciSerialization` AST evaluator (`eval_coeff`, src/rule_eval.jl).

Declarative companion: `discretizations/grids/duo/faq/primal_geometry.esm` (D2a /
esd-heg.6). That document expresses the per-cell geometry of `build_duo_grid`
(src/grids/duo.jl) — the cell-center centroid (normalized mean of the three corner
unit vectors), the geographic `lon`/`lat` (`atan2`/`asin`), and the spherical-excess
cell area (L'Huilier) — as the RFC `semiring-faq-unified-ir` §8.1 geometry-FAQ
pattern: per cell, a `sum_product` aggregate over the three corner vertices producing
each scalar metric from the inverse-trig / `sqrt` / `tan` scalar leaves.

This file is the thin **evaluation bridge**: it shapes the unit vertices + faces into
the per-cell scalar inputs and routes every piece of grid ARITHMETIC through ESS's
single AST pathway (`eval_coeff`). ESD hosts no shadow evaluator (AGENTS.md
single-pathway invariant), so the determinism contract — and the cross-binding
byte-identity of the output — lives entirely in `EarthSciSerialization`. The only
ESD-side logic is gathering the corner coordinates and the `Float64(...)` / `T(...)`
coercions; the float arithmetic is all the rule AST.

The transcription mirrors `build_duo_grid`'s per-cell loop bit-for-bit (T=Float64):
squares are products (`x*x`, not `^`); the centroid sum, the squared norm, and the
three side-arc dot products fold in corner / space order; the L'Huilier tan-product
is left-folded; `area = (4·atan(√t)) · (R·R)` keeps the imperative association
`a_unit * R2`. The single intentional divergence is dropping the `clamp` / `max`
clipping guards: on a valid icosahedral mesh every `acos` / `asin` argument and the
L'Huilier radicand already lie in range, so the guards never fire and the bytes match
(proven 0-ULP in test/test_duo_primal_geometry_faq.jl across the level ladder).
"""

# AST nodes mirroring the equation bodies of
# `discretizations/grids/duo/faq/primal_geometry.esm`, lifted verbatim from
# test/test_duo_primal_geometry_faq.jl. Held as module consts so per-cell evaluation
# does not rebuild the trees; `eval_coeff` (src/rule_eval.jl) re-parses each call —
# the same single pathway the rule catalog uses. Scalar corner leaves: a1..a3 (corner
# 1 = face vertex 1), b1..b3 (corner 2), c1..c3 (corner 3), and the radius R.
_pg_mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

# centroid: m = a+b+c per component; n = sqrt(mx*mx+my*my+mz*mz); u = m/n.
const _PG_MX = _pg_mk("+", "a1", "b1", "c1")
const _PG_MY = _pg_mk("+", "a2", "b2", "c2")
const _PG_MZ = _pg_mk("+", "a3", "b3", "c3")
const _PG_NSQ = _pg_mk("+",
    _pg_mk("*", _PG_MX, _PG_MX), _pg_mk("*", _PG_MY, _PG_MY), _pg_mk("*", _PG_MZ, _PG_MZ))   # squares as *, not ^
const _PG_NRM = _pg_mk("sqrt", _PG_NSQ)
const _PG_UX = _pg_mk("/", _PG_MX, _PG_NRM)
const _PG_UY = _pg_mk("/", _PG_MY, _PG_NRM)
const _PG_UZ = _pg_mk("/", _PG_MZ, _PG_NRM)

# cell-center cartesian (R-scaled), geographic lon/lat of the unit centroid.
const _PG_CART_X = _pg_mk("*", "R", _PG_UX)
const _PG_CART_Y = _pg_mk("*", "R", _PG_UY)
const _PG_CART_Z = _pg_mk("*", "R", _PG_UZ)
const _PG_LON = _pg_mk("atan2", _PG_UY, _PG_UX)
const _PG_LAT = _pg_mk("asin", _PG_UZ)                                    # no clamp

# area (L'Huilier): da=acos(b·c), db=acos(c·a), dc=acos(a·b); s=(da+db+dc)/2;
# t=tan(s/2)·tan((s-da)/2)·tan((s-db)/2)·tan((s-dc)/2); area = (4·atan(√t)) · (R·R).
const _PG_DOT_BC = _pg_mk("+", _pg_mk("*", "b1", "c1"), _pg_mk("*", "b2", "c2"), _pg_mk("*", "b3", "c3"))
const _PG_DOT_CA = _pg_mk("+", _pg_mk("*", "c1", "a1"), _pg_mk("*", "c2", "a2"), _pg_mk("*", "c3", "a3"))
const _PG_DOT_AB = _pg_mk("+", _pg_mk("*", "a1", "b1"), _pg_mk("*", "a2", "b2"), _pg_mk("*", "a3", "b3"))
const _PG_DA = _pg_mk("acos", _PG_DOT_BC)
const _PG_DB = _pg_mk("acos", _PG_DOT_CA)
const _PG_DC = _pg_mk("acos", _PG_DOT_AB)
const _PG_S = _pg_mk("*", 0.5, _pg_mk("+", _PG_DA, _PG_DB, _PG_DC))
_pg_half(x) = _pg_mk("/", x, 2.0)
const _PG_T = _pg_mk("*",
    _pg_mk("tan", _pg_half(_PG_S)),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DA))),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DB))),
    _pg_mk("tan", _pg_half(_pg_mk("-", _PG_S, _PG_DC))))
const _PG_AREA = _pg_mk("*", _pg_mk("*", 4.0, _pg_mk("atan", _pg_mk("sqrt", _PG_T))), _pg_mk("*", "R", "R"))   # no max guard

"""
    primal_geometry_faq(::Type{T}, vertices_unit, faces, R) where {T}
        -> (cell_cart, lon, lat, area)

Materialize the DUO primal-cell geometry from the declarative FAQ
(`discretizations/grids/duo/faq/primal_geometry.esm`), routing all arithmetic through
the landed ESS AST evaluator (`eval_coeff`).

Arguments:
- `vertices_unit::AbstractMatrix{<:Real}` — `(3, Nv)` cartesian **unit-sphere** vertex
  coordinates (NOT R-scaled), as produced by `_subdivide_icosahedron`.
- `faces::AbstractMatrix{<:Integer}` — `(3, Nc)` 1-based vertex indices per triangular
  cell.
- `R::Real` — sphere radius.

Returns a `NamedTuple`:
- `cell_cart :: Matrix{T}` `(3, Nc)` — R-scaled cell-center cartesian coordinates.
- `lon :: Vector{T}` `(Nc,)` — cell-center longitude (rad).
- `lat :: Vector{T}` `(Nc,)` — cell-center latitude (rad).
- `area :: Vector{T}` `(Nc,)` — spherical-triangle area (R² scaled).

Byte-identical (T=Float64) to `build_duo_grid`'s cell-geometry loop (src/grids/duo.jl):
the per-cell bindings and ASTs mirror the imperative float ops exactly, and the
`clamp` / `max` guards are dropped because they never fire on a valid mesh.
"""
function primal_geometry_faq(::Type{T}, vertices_unit::AbstractMatrix{<:Real},
        faces::AbstractMatrix{<:Integer}, R::Real) where {T}
    V = vertices_unit
    F = faces
    Nc = size(F, 2)
    Rf = Float64(R)

    cell_cart = Matrix{T}(undef, 3, Nc)
    lon = Vector{T}(undef, Nc)
    lat = Vector{T}(undef, Nc)
    area = Vector{T}(undef, Nc)

    @inbounds for c in 1:Nc
        a_i = F[1, c]; b_i = F[2, c]; c_i = F[3, c]
        # One bindings Dict per cell drives every co-derived metric (centroid →
        # cart/lon/lat, and the L'Huilier area), mirroring the single imperative loop.
        b = Dict{String, Float64}(
            "a1" => Float64(V[1, a_i]), "a2" => Float64(V[2, a_i]), "a3" => Float64(V[3, a_i]),
            "b1" => Float64(V[1, b_i]), "b2" => Float64(V[2, b_i]), "b3" => Float64(V[3, b_i]),
            "c1" => Float64(V[1, c_i]), "c2" => Float64(V[2, c_i]), "c3" => Float64(V[3, c_i]),
            "R" => Rf,
        )
        cell_cart[1, c] = T(eval_coeff(_PG_CART_X, b))
        cell_cart[2, c] = T(eval_coeff(_PG_CART_Y, b))
        cell_cart[3, c] = T(eval_coeff(_PG_CART_Z, b))
        lon[c] = T(eval_coeff(_PG_LON, b))
        lat[c] = T(eval_coeff(_PG_LAT, b))
        area[c] = T(eval_coeff(_PG_AREA, b))
    end

    return (cell_cart = cell_cart, lon = lon, lat = lat, area = area)
end
