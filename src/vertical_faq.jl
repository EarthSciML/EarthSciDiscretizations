"""
    vertical_faq.jl — vertical (1-D column) structured-grid construction via the
    landed M1 elementwise FAQ path (`EarthSciSerialization` AST evaluator).

Declarative companion: `discretizations/grids/vertical/rules/vertical_construction.esm`
expresses this construction as a semiring-FAQ document (RFC `semiring-faq-unified-ir`
§5.1/§5.2: two interval index-sets — `interfaces` (nz+1 half-levels) and `layers`
(nz cells) — plus elementwise `arrayop` bodies; the M1 path, NO value-invention, all
CONST cadence). It reuses the structured-grid FAQ template established by
`src/cartesian_faq.jl` (`cartesian_construction.esm`, esd-3we.1 / S1).

This module is the thin **evaluation bridge**: it routes every piece of vertical
grid ARITHMETIC — the per-coordinate level synthesis, the cell-center/width
derivations, the cell-volume, and the named layer metrics — through ESS's AST
evaluator (`eval_coeff`, the single pathway in `src/rule_eval.jl`). ESD hosts no
shadow evaluator (AGENTS.md single-pathway invariant; GRIDS_API §4.3), so the
determinism contract — and therefore the cross-binding byte-identity of the output —
lives entirely in `EarthSciSerialization` (`CONFORMANCE_SPEC.md` §5.5). The only
ESD-side logic is the *structural* part — the 0-sentinel neighbor maps and the
boundary masks — exactly as `cartesian_faq.jl` materializes its neighbor / boundary
arrays structurally while routing the affine/elementwise arithmetic through ESS.

Level synthesis mirrors the imperative `_vertical` generator (`src/grids/vertical.jl`)
branch-for-branch so the FAQ output is bit-identical to its trait arrays:

- `:eta` — synthesized hybrid sigma-pressure `levels[k] = ak[k]/p0 + bk[k]` (ak/bk
  are DATA inputs); one divide-add per interface, routed through ESS.
- `:sigma` / `:hybrid_sigma_theta` — the uniform-sigma affine map `levels[k] =
  1 - (k-1)/nz` routed through ESS when the grid was built uniform; a grid built
  from explicit (non-uniform) sigma `levels` instead carries them as DATA (the
  supplied-array branch, the analogue of cartesian's non-uniform supplied `edges`).
- `:z` / `:theta` / `:z_star` — interface `levels` supplied verbatim as DATA.

The downstream derivations (centers, widths, volume, metrics) are pure arithmetic on
those levels and always ride ESS.
"""

# AST nodes mirroring the equation bodies of
# `discretizations/grids/vertical/rules/vertical_construction.esm`. Held as module
# consts so the per-element evaluation does not re-allocate the trees; `eval_coeff`
# (src/rule_eval.jl) re-parses each call, the same single pathway the rule catalog uses.

# (A) uniform sigma level: levels[k] = 1 - (k-1)/nz  (binding `k` is the 0-based
# interface index, i.e. (1-based i) - 1; matches `one(T) - T(k)/T(nz)`).
const _VERT_SIGMA_LEVEL_NODE = Dict{String, Any}(
    "op" => "-",
    "args" => Any[1, Dict{String, Any}("op" => "/", "args" => Any["k", "nz"])],
)
# (A) eta hybrid level: levels[k] = ak[k]/p0 + bk[k]  (matches `ak ./ p0 .+ bk`).
const _VERT_ETA_LEVEL_NODE = Dict{String, Any}(
    "op" => "+",
    "args" => Any[Dict{String, Any}("op" => "/", "args" => Any["ak", "p0"]), "bk"],
)
# (B/C) two-point average: (a + b)/2 — the midpoint reduction. Used for cell
# `centers` ((levels[k]+levels[k+1])/2) and the metric layer-averages
# (metric_ak/metric_bk = (coeff[k]+coeff[k+1])/2).
const _VERT_AVG2_NODE = Dict{String, Any}(
    "op" => "/",
    "args" => Any[Dict{String, Any}("op" => "+", "args" => Any["a", "b"]), 2],
)
# (B) layer width: abs(b - a) = abs(levels[k+1] - levels[k]) — positive thickness
# for decreasing (sigma) and increasing (z/theta) levels alike.
const _VERT_WIDTH_NODE = Dict{String, Any}(
    "op" => "abs",
    "args" => Any[Dict{String, Any}("op" => "-", "args" => Any["b", "a"])],
)
# (C) layer-average reference pressure: (p_lo + p_hi)/2 with
# p_lo = ak_lo + bk_lo*p0 and p_hi = ak_hi + bk_hi*p0 — matches the imperative
# `(g.ak[k] + g.bk[k]*g.p0 + g.ak[k+1] + g.bk[k+1]*g.p0)/2` evaluation order.
const _VERT_PRESSURE_NODE = Dict{String, Any}(
    "op" => "/",
    "args" => Any[
        Dict{String, Any}(
            "op" => "+",
            "args" => Any[
                Dict{String, Any}(
                    "op" => "+",
                    "args" => Any["ak_lo", Dict{String, Any}("op" => "*", "args" => Any["bk_lo", "p0"])]
                ),
                Dict{String, Any}(
                    "op" => "+",
                    "args" => Any["ak_hi", Dict{String, Any}("op" => "*", "args" => Any["bk_hi", "p0"])]
                ),
            ],
        ),
        2,
    ],
)

"Bit-exact equality of two Float64-coercible vectors (distinguishes ±0.0 / NaN payloads)."
function _vert_bitequal(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        reinterpret(UInt64, Float64(a[i])) == reinterpret(UInt64, Float64(b[i])) || return false
    end
    return true
end

"""
    _faq_levels(grid::VerticalGrid{T}) -> Vector{T}

Synthesize the interface (half-level) coordinates via the M1 FAQ, branching on the
coordinate kind exactly as the imperative `_vertical` generator does. Eta levels are
the divide-add `ak[k]/p0 + bk[k]` routed through ESS; uniform sigma is the affine map
`1 - (k-1)/nz` routed through ESS (and used only when it reproduces the grid's stored
levels bit-for-bit — a grid built from explicit sigma levels falls back to that
supplied DATA); every other kind supplies its `levels` verbatim. The result equals
`grid.levels` bit-for-bit, so the bridge stays anchored to the imperative builder.
"""
function _faq_levels(grid::VerticalGrid{T}) where {T}
    coord = grid.coordinate
    nz = n_cells(grid)
    if coord === :eta
        p0 = Float64(grid.p0)
        return T[
            T(
                    eval_coeff(
                        _VERT_ETA_LEVEL_NODE,
                        Dict("ak" => Float64(grid.ak[k]), "p0" => p0, "bk" => Float64(grid.bk[k]))
                    )
                )
                for k in 1:(nz + 1)
        ]
    elseif coord === :sigma || coord === :hybrid_sigma_theta
        cand = T[
            T(
                    eval_coeff(
                        _VERT_SIGMA_LEVEL_NODE,
                        Dict("k" => Float64(k - 1), "nz" => Float64(nz))
                    )
                )
                for k in 1:(nz + 1)
        ]
        return _vert_bitequal(cand, grid.levels) ? cand : copy(grid.levels)
    else  # :z, :theta, :z_star — supplied verbatim as DATA
        return copy(grid.levels)
    end
end

"""
    _faq_centers_widths(levels::Vector{T}) -> (centers, widths)

Elementwise midpoint / absolute-difference over the `layers` interval, routed through
ESS: `centers[k] = (levels[k] + levels[k+1])/2`, `widths[k] = abs(levels[k+1] -
levels[k])`.
"""
function _faq_centers_widths(levels::Vector{T}) where {T}
    n = length(levels) - 1
    centers = Vector{T}(undef, n)
    widths = Vector{T}(undef, n)
    @inbounds for k in 1:n
        b = Dict("a" => Float64(levels[k]), "b" => Float64(levels[k + 1]))
        centers[k] = T(eval_coeff(_VERT_AVG2_NODE, b))
        widths[k] = T(eval_coeff(_VERT_WIDTH_NODE, b))
    end
    return centers, widths
end

"""
    vertical_construction_faq(grid::VerticalGrid{T}) -> NamedTuple

Materialize a vertical column's full construction — interface levels, cell centers /
widths / volume, the named layer metrics, the plus/minus-k neighbor maps and the
boundary masks — from the declarative M1 elementwise FAQ, evaluating all arithmetic
through the landed ESS AST evaluator. Returns a NamedTuple whose fields match the
imperative `VerticalGrid` trait arrays (`src/grids/vertical.jl`) so the conformance
harness can assert byte/ULP equality:

- `levels` :: `Vector{T}` — interface coordinates (length nz+1), synthesized per
  coordinate kind (`_faq_levels`); equals `half_levels(grid)`.
- `centers`, `widths` :: `Vector{T}` — cell centers / thicknesses (length nz),
  matching `cell_centers(grid)` / `cell_widths(grid)` (and `layer_thickness`).
- `cell_volume` :: `Vector{T}` — per-layer measure (= `widths`), matching
  `cell_volume(grid)`.
- `metric_dz`, `metric_z` :: `Vector{T}` — `metric_eval(grid, :dz|:z, ·)` over the
  column. `metric_sigma`, `metric_pressure`, `metric_ak`, `metric_bk` :: either a
  `Vector{T}` (when the metric is defined for this coordinate) or `nothing`, matching
  the imperative `metric_eval` availability (sigma for sigma-like kinds; pressure/ak/bk
  for grids carrying hybrid coefficients).
- `neighbor_minus`, `neighbor_plus` :: `Vector{Int}` — the offset −1 / +1 layer
  neighbor with the `0` sentinel for off-column, matching `neighbor_indices(grid, :z, ∓1)`.
- `boundary_lower`, `boundary_upper` :: `Vector{Bool}` — bottom / top layer, matching
  `boundary_mask(grid, :z, :lower/:upper)`.
- `half_levels`, `layer_thickness` :: aliases of `levels` / `widths`; `ak`, `bk`, `p0`
  — the `pressure_coefficients(grid)` (ak/bk empty when the kind has none).

The level/center/width/metric arithmetic rides ESS's determinism contract; the
structural integer arrays (neighbor linearization, masks) are pure index logic, the
1-D analogue of the structural materialization in `cartesian_construction_faq`.
"""
function vertical_construction_faq(grid::VerticalGrid{T}) where {T}
    levels = _faq_levels(grid)
    centers, widths = _faq_centers_widths(levels)
    n = length(widths)

    # cell_volume / metric :dz / :z / :sigma are identity reads of the already
    # ESS-derived widths / centers (matching the imperative accessors).
    cell_volume = copy(widths)
    metric_dz = copy(widths)
    metric_z = copy(centers)

    sigma_like = grid.coordinate === :sigma ||
        grid.coordinate === :hybrid_sigma_theta ||
        grid.coordinate === :eta
    metric_sigma = sigma_like ? copy(centers) : nothing

    has_ak = length(grid.ak) > 0
    has_bk = length(grid.bk) > 0
    p0 = Float64(grid.p0)

    metric_pressure = if has_ak && has_bk
        T[
            T(
                    eval_coeff(
                        _VERT_PRESSURE_NODE, Dict(
                            "ak_lo" => Float64(grid.ak[k]), "bk_lo" => Float64(grid.bk[k]),
                            "ak_hi" => Float64(grid.ak[k + 1]), "bk_hi" => Float64(grid.bk[k + 1]),
                            "p0" => p0
                        )
                    )
                ) for k in 1:n
        ]
    else
        nothing
    end
    metric_ak = has_ak ?
        T[
            T(
                eval_coeff(
                    _VERT_AVG2_NODE,
                    Dict("a" => Float64(grid.ak[k]), "b" => Float64(grid.ak[k + 1]))
                )
            ) for k in 1:n
        ] :
        nothing
    metric_bk = has_bk ?
        T[
            T(
                eval_coeff(
                    _VERT_AVG2_NODE,
                    Dict("a" => Float64(grid.bk[k]), "b" => Float64(grid.bk[k + 1]))
                )
            ) for k in 1:n
        ] :
        nothing

    # --- structural arrays (0-sentinel neighbor maps, boundary masks) -----------
    neighbor_minus = Int[k > 1 ? k - 1 : 0 for k in 1:n]
    neighbor_plus = Int[k < n ? k + 1 : 0 for k in 1:n]
    boundary_lower = Bool[k == 1 for k in 1:n]
    boundary_upper = Bool[k == n for k in 1:n]

    return (
        levels = levels,
        centers = centers,
        widths = widths,
        cell_volume = cell_volume,
        metric_dz = metric_dz,
        metric_z = metric_z,
        metric_sigma = metric_sigma,
        metric_pressure = metric_pressure,
        metric_ak = metric_ak,
        metric_bk = metric_bk,
        neighbor_minus = neighbor_minus,
        neighbor_plus = neighbor_plus,
        boundary_lower = boundary_lower,
        boundary_upper = boundary_upper,
        half_levels = levels,
        layer_thickness = widths,
        ak = copy(grid.ak),
        bk = copy(grid.bk),
        p0 = grid.p0,
    )
end
