"""
Vertical grid family (1D column).

Implements the public API contract defined in `docs/GRIDS_API.md` §2.3, §3.1
and §7 for the Julia binding. The generator entry point lives in the
`EarthSciDiscretizations.grids` submodule (wired at the bottom of
`src/EarthSciDiscretizations.jl`) and is invoked as:

    using EarthSciDiscretizations
    g = EarthSciDiscretizations.grids.vertical(; coordinate = :sigma, nz = 16)
    g = EarthSciDiscretizations.grids.vertical(; coordinate = :z,
                                                 levels = [0.0, 500.0, 1500.0])
    g = EarthSciDiscretizations.grids.vertical(;
            coordinate = :eta, ak = ak_vec, bk = bk_vec, p0 = 1.0e5)

Supported coordinate kinds (rhymes with the Python sibling,
`python/src/earthsci_toolkit/grids/vertical.py`):

- `:sigma` — terrain-following, dimensionless ∈ [0, 1]; 1 at surface, 0 at top.
- `:eta` — hybrid sigma-pressure (NCAR / CAM style). Synthesized sigma
  = `ak/p0 + bk` at interfaces; must be strictly decreasing.
- `:z` — geometric altitude (m); strictly increasing.
- `:theta` — potential temperature (K); strictly increasing.
- `:hybrid_sigma_theta` — blended sigma → theta; sigma-like levels.
- `:z_star` — generalized geometric height; strictly increasing.

Per the 2026-04-20 scope correction on the vertical fixtures bead, the
`.esm` lowering is a small declarative config (family + coordinate kind +
interface levels + hybrid coefficients when applicable + provenance), NOT a
serialized geometry blob. Cell centers and widths are derived on demand
from `levels` via pure arithmetic.
"""

# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

struct VerticalGrid{T <: AbstractFloat} <: AbstractGrid
    coordinate::Symbol
    levels::Vector{T}
    centers::Vector{T}
    widths::Vector{T}
    ak::Vector{T}
    bk::Vector{T}
    p0::T
    dtype::String
    ghosts::Int
end

Base.ndims(::VerticalGrid) = 1
Base.eltype(::VerticalGrid{T}) where {T} = T

n_cells(g::VerticalGrid) = length(g.levels) - 1
n_vertices(g::VerticalGrid) = length(g.levels)
n_edges(g::VerticalGrid) = length(g.levels) - 1

const _VERTICAL_FAMILY_VERSION = "1.0.0"

const _VERTICAL_COORDINATES =
    (:sigma, :eta, :z, :theta, :hybrid_sigma_theta, :z_star)

const _VERTICAL_METRIC_NAMES = (:dz, :z, :sigma, :pressure, :ak, :bk)

# ---------------------------------------------------------------------------
# Parse + validation helpers
# ---------------------------------------------------------------------------

_parse_vertical_coordinate(c::Symbol) =
    c in _VERTICAL_COORDINATES ? c :
        throw(ArgumentError(
            "vertical: unknown coordinate :$c; expected one of " *
                "$(_VERTICAL_COORDINATES)"
        ))
_parse_vertical_coordinate(c::AbstractString) =
    _parse_vertical_coordinate(Symbol(c))

function _vertical_uniform_sigma(nz::Int, ::Type{T}) where {T}
    nz ≥ 1 || throw(ArgumentError("vertical: nz must be ≥ 1; got $nz"))
    return T[one(T) - T(k) / T(nz) for k in 0:nz]
end

function _vertical_coerce_levels(
        levels, ::Type{T};
        must_decrease::Bool,
        domain::Union{Tuple{Real, Real}, Nothing} = nothing,
        label::AbstractString = "levels"
    ) where {T}
    arr = T[T(x) for x in levels]
    length(arr) ≥ 2 ||
        throw(ArgumentError("vertical: `$label` must have ≥ 2 entries; got $(length(arr))"))
    for x in arr
        isfinite(x) ||
            throw(DomainError(x, "vertical: `$label` must be finite"))
    end
    if domain !== nothing
        lo, hi = T(domain[1]), T(domain[2])
        for x in arr
            (x < lo || x > hi) &&
                throw(DomainError(x,
                    "vertical: `$label` entries must lie in [$lo, $hi]"))
        end
    end
    if must_decrease
        for k in 1:(length(arr) - 1)
            arr[k + 1] < arr[k] ||
                throw(DomainError(arr,
                    "vertical: `$label` must be strictly decreasing"))
        end
    else
        for k in 1:(length(arr) - 1)
            arr[k + 1] > arr[k] ||
                throw(DomainError(arr,
                    "vertical: `$label` must be strictly increasing"))
        end
    end
    return arr
end

function _vertical_coerce_hybrid(
        coeffs, ::Type{T}, expected_len::Int, label::AbstractString
    ) where {T}
    coeffs === nothing &&
        throw(ArgumentError("vertical: `$label` is required"))
    arr = T[T(x) for x in coeffs]
    length(arr) == expected_len ||
        throw(ArgumentError(
            "vertical: `$label` must have length nz+1 = $expected_len; " *
                "got $(length(arr))"
        ))
    for x in arr
        isfinite(x) ||
            throw(DomainError(x, "vertical: `$label` must be finite"))
    end
    return arr
end

function _vertical_centers_and_widths(levels::Vector{T}) where {T}
    n = length(levels) - 1
    centers = Vector{T}(undef, n)
    widths = Vector{T}(undef, n)
    for k in 1:n
        centers[k] = (levels[k] + levels[k + 1]) / T(2)
        widths[k] = abs(levels[k + 1] - levels[k])
    end
    return centers, widths
end

# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

"""
    _vertical(; coordinate, nz=nothing, levels=nothing,
                ak=nothing, bk=nothing, p0=1.0e5, transition=nothing,
                dtype=Float64, ghosts=0) -> VerticalGrid{T}

Primary constructor for the vertical family. See `docs/GRIDS_API.md` §2.3 for
the cross-binding signature contract. Public entry point (Julia):
`EarthSciDiscretizations.grids.vertical`.
"""
function _vertical(;
        coordinate = nothing,
        nz::Union{Int, Nothing} = nothing,
        levels = nothing,
        ak = nothing,
        bk = nothing,
        p0::Real = 1.0e5,
        transition::Union{Real, Nothing} = nothing,
        dtype::Type{T} = Float64,
        ghosts::Int = 0
    ) where {T <: AbstractFloat}

    coordinate === nothing &&
        throw(ArgumentError("vertical: required option `coordinate` is missing"))
    coord = _parse_vertical_coordinate(coordinate)

    ghosts ≥ 0 || throw(ArgumentError("vertical: ghosts must be ≥ 0; got $ghosts"))

    if nz !== nothing
        nz ≥ 1 || throw(ArgumentError("vertical: nz must be ≥ 1; got $nz"))
    end

    p0_T = T(p0)
    (p0_T > zero(T) && isfinite(p0_T)) ||
        throw(DomainError(p0, "vertical: p0 must be positive and finite"))

    local lv::Vector{T}
    local ak_arr::Vector{T}
    local bk_arr::Vector{T}

    if coord === :sigma
        if levels !== nothing
            lv = _vertical_coerce_levels(
                levels, T; must_decrease = true, domain = (0.0, 1.0)
            )
            nz !== nothing && length(lv) != nz + 1 &&
                throw(ArgumentError(
                    "vertical: sigma: nz=$nz inconsistent with levels length $(length(lv))"
                ))
        else
            nz === nothing &&
                throw(ArgumentError("vertical: sigma requires `nz` or `levels`"))
            lv = _vertical_uniform_sigma(nz, T)
        end
        ak_arr = T[]
        bk_arr = T[]

    elseif coord === :z || coord === :theta || coord === :z_star
        levels === nothing &&
            throw(ArgumentError("vertical: :$coord requires explicit `levels`"))
        lv = _vertical_coerce_levels(levels, T; must_decrease = false)
        nz !== nothing && length(lv) != nz + 1 &&
            throw(ArgumentError(
                "vertical: $coord: nz=$nz inconsistent with levels length $(length(lv))"
            ))
        ak_arr = T[]
        bk_arr = T[]

    elseif coord === :eta
        ak === nothing && throw(ArgumentError("vertical: :eta requires `ak` (length nz+1)"))
        bk === nothing && throw(ArgumentError("vertical: :eta requires `bk` (length nz+1)"))
        ak_probe = T[T(x) for x in ak]
        bk_probe = T[T(x) for x in bk]
        length(ak_probe) == length(bk_probe) ||
            throw(ArgumentError(
                "vertical: :eta ak/bk must have equal length; " *
                    "got $(length(ak_probe)) vs $(length(bk_probe))"
            ))
        nz_eff = nz === nothing ? length(ak_probe) - 1 : nz
        nz_eff ≥ 1 || throw(ArgumentError("vertical: nz must be ≥ 1; got $nz_eff"))
        ak_arr = _vertical_coerce_hybrid(ak, T, nz_eff + 1, "ak")
        bk_arr = _vertical_coerce_hybrid(bk, T, nz_eff + 1, "bk")
        sigma = ak_arr ./ p0_T .+ bk_arr
        for k in 1:(length(sigma) - 1)
            sigma[k + 1] < sigma[k] ||
                throw(DomainError(sigma,
                    "vertical: :eta synthesized sigma (ak/p0 + bk) must be " *
                        "strictly decreasing"
                ))
        end
        lv = sigma

    elseif coord === :hybrid_sigma_theta
        if levels === nothing && nz === nothing
            throw(ArgumentError(
                "vertical: :hybrid_sigma_theta requires `nz` or `levels`"
            ))
        end
        if levels !== nothing
            lv = _vertical_coerce_levels(
                levels, T; must_decrease = true, domain = (0.0, 1.0)
            )
            nz !== nothing && length(lv) != nz + 1 &&
                throw(ArgumentError(
                    "vertical: hybrid_sigma_theta: nz=$nz inconsistent with " *
                        "levels length $(length(lv))"
                ))
        else
            lv = _vertical_uniform_sigma(nz, T)
        end
        if transition !== nothing
            (zero(T) < T(transition) < one(T)) ||
                throw(DomainError(transition,
                    "vertical: hybrid_sigma_theta `transition` must be in (0, 1)"
                ))
        end
        ak_arr = ak === nothing ? T[] :
            _vertical_coerce_hybrid(ak, T, length(lv), "ak")
        bk_arr = bk === nothing ? T[] :
            _vertical_coerce_hybrid(bk, T, length(lv), "bk")

    else
        # Unreachable: _parse_vertical_coordinate already validated the set.
        throw(AssertionError("vertical: unreachable coordinate :$coord"))
    end

    centers, widths = _vertical_centers_and_widths(lv)
    return VerticalGrid{T}(
        coord, lv, centers, widths, ak_arr, bk_arr, p0_T,
        _dtype_string(T), ghosts,
    )
end

# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------

function _vertical_check_layer(g::VerticalGrid, k::Int)
    (1 ≤ k ≤ n_cells(g)) ||
        throw(BoundsError(g, k))
    return nothing
end

"""
    cell_centers(g::VerticalGrid, [k::Int]) -> Vector{T} | T

With no index, return all mid-layer coordinate values (length `nz`). With an
integer `k` (1-based), return the scalar center at that layer.
"""
cell_centers(g::VerticalGrid) = g.centers
function cell_centers(g::VerticalGrid, k::Integer)
    _vertical_check_layer(g, Int(k))
    return g.centers[k]
end

"""
    cell_widths(g::VerticalGrid, [k::Int]) -> Vector{T} | T

Layer thicknesses in native units, always positive. With no index returns
the full widths vector (length `nz`).
"""
cell_widths(g::VerticalGrid) = g.widths
function cell_widths(g::VerticalGrid, k::Integer)
    _vertical_check_layer(g, Int(k))
    return g.widths[k]
end

"""
    neighbors(g::VerticalGrid, k::Int) -> Dict{Symbol,Int}

Axis-aligned vertical neighbours of layer `k` (1-based). Returns a dict
with `:down` → `k-1` and/or `:up` → `k+1`. Bottom / top layers drop the
out-of-range side, keeping the return type concrete.
"""
function neighbors(g::VerticalGrid, k::Integer)
    _vertical_check_layer(g, Int(k))
    out = Dict{Symbol, Int}()
    if k > 1
        out[:down] = Int(k) - 1
    end
    if Int(k) < n_cells(g)
        out[:up] = Int(k) + 1
    end
    return out
end

"""
    metric_eval(g::VerticalGrid, name::Symbol, k::Int) -> T

Evaluate a named metric at layer `k` (1-based). Supported names:

- `:dz` — layer thickness (native units).
- `:z` — cell-centre value (native units).
- `:sigma` — sigma at cell centre; valid for sigma-like coordinates
  (`:sigma`, `:hybrid_sigma_theta`, `:eta`).
- `:pressure` — reference pressure at cell centre, averaged across the
  layer's two interfaces (`p = ak + bk * p0`). Requires hybrid
  coefficients.
- `:ak`, `:bk` — hybrid coefficients averaged across the layer's two
  interfaces. Requires hybrid coefficients.
"""
function metric_eval(g::VerticalGrid{T}, name::Symbol, k::Integer) where {T}
    name in _VERTICAL_METRIC_NAMES ||
        throw(ArgumentError("vertical: unknown metric :$name"))
    _vertical_check_layer(g, Int(k))
    if name === :dz
        return g.widths[k]
    elseif name === :z
        return g.centers[k]
    elseif name === :sigma
        (g.coordinate === :sigma ||
            g.coordinate === :hybrid_sigma_theta ||
            g.coordinate === :eta) ||
            throw(ArgumentError(
                "vertical: :sigma undefined for coordinate :$(g.coordinate)"
            ))
        return g.centers[k]
    elseif name === :pressure
        (length(g.ak) > 0 && length(g.bk) > 0) ||
            throw(ArgumentError(
                "vertical: :pressure requires hybrid ak/bk " *
                    "(coordinate :$(g.coordinate) has none)"
            ))
        p_lo = g.ak[k] + g.bk[k] * g.p0
        p_hi = g.ak[k + 1] + g.bk[k + 1] * g.p0
        return (p_lo + p_hi) / T(2)
    elseif name === :ak
        length(g.ak) > 0 ||
            throw(ArgumentError("vertical: :ak unavailable (no hybrid coefficients)"))
        return (g.ak[k] + g.ak[k + 1]) / T(2)
    else  # :bk
        length(g.bk) > 0 ||
            throw(ArgumentError("vertical: :bk unavailable (no hybrid coefficients)"))
        return (g.bk[k] + g.bk[k + 1]) / T(2)
    end
end

metric_eval(g::VerticalGrid, name::AbstractString, k::Integer) =
    metric_eval(g, Symbol(name), k)

# ---------------------------------------------------------------------------
# .esm lowering + provenance
# ---------------------------------------------------------------------------

function _vertical_provenance(g::VerticalGrid)
    return Dict{String, Any}(
        "binding" => "julia",
        "binding_version" => _VERTICAL_FAMILY_VERSION,
        "family" => "vertical",
        "version" => _VERTICAL_FAMILY_VERSION,
        "coordinate" => String(g.coordinate),
        "dtype" => g.dtype,
    )
end

"""
    to_esm(g::VerticalGrid) -> Dict{String,Any}

Declarative `.esm` lowering per the 2026-04-20 scope correction. The output
carries the family + coordinate kind + interface levels (+ hybrid
coefficients when applicable) + provenance; derived arrays (centres,
widths) do not appear in the wire form. Matches the Python reference
binding's key set so cross-binding byte-equality is achievable by the
conformance harness.
"""
function to_esm(g::VerticalGrid)
    options = Dict{String, Any}(
        "coordinate" => String(g.coordinate),
        "nz" => n_cells(g),
        "levels" => Float64[Float64(x) for x in g.levels],
    )
    if length(g.ak) > 0
        options["ak"] = Float64[Float64(x) for x in g.ak]
    end
    if length(g.bk) > 0
        options["bk"] = Float64[Float64(x) for x in g.bk]
    end
    if g.coordinate === :eta || g.coordinate === :hybrid_sigma_theta ||
            length(g.ak) > 0 || length(g.bk) > 0
        options["p0"] = Float64(g.p0)
    end
    return Dict{String, Any}(
        "family" => "vertical",
        "topology" => "column",
        "dtype" => g.dtype,
        "ndim" => 1,
        "ghosts" => g.ghosts,
        "n_cells" => n_cells(g),
        "n_vertices" => n_vertices(g),
        "n_edges" => n_edges(g),
        "options" => options,
        "provenance" => _vertical_provenance(g),
        "schema_version" => _VERTICAL_FAMILY_VERSION,
    )
end

family(::VerticalGrid) = "vertical"

# Public namespace registration: the `vertical` entry point in the
# `EarthSciDiscretizations.grids` submodule (defined at the bottom of
# `src/EarthSciDiscretizations.jl`) aliases `_vertical` here. See
# GRIDS_API.md §2.3 for the public call form.
