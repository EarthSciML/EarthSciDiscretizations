"""
Lat-lon grid family (regular + reduced-Gaussian variants).

Implements the public API contract defined in `docs/GRIDS_API.md` §2.3 / §3.1 /
§7 for the Julia binding. The generator entry point lives in the
`EarthSciDiscretizations.grids` submodule (wired at the bottom of
`src/EarthSciDiscretizations.jl`) and is invoked as:

    using EarthSciDiscretizations
    g = EarthSciDiscretizations.grids.lat_lon(; nlon = 72, nlat = 36)
    g = EarthSciDiscretizations.grids.lat_lon(;
            variant = :reduced_gaussian,
            nlon_per_row = [4, 8, 12, 12, 8, 4],
            R = 1.0,
        )

Two variants are supported:

- `:regular`          — strictly uniform in lon and lat.
- `:reduced_gaussian` — per-row `nlon_per_row` with optional user-supplied
  `lat_edges` for genuine Gaussian-quadrature latitudes. When `lat_edges` is
  omitted it defaults to equal-angle edges, which is a test-only convenience.

Polar-singularity handling is declared via `pole_policy` but only `:none` is
implemented in this phase; `:average` / `:fold` raise `ArgumentError` at
construction so downstream callers see the not-yet-implemented surface rather
than silently incorrect geometry.

Per the 2026-04-20 scope correction, the `.esm` lowering is a small declarative
config (family + variant + generator + params + provenance), NOT a serialized
geometry blob. Cell-center coords, neighbours, and metrics are derived on
demand by the accessors below.
"""

const _LATLON_FAMILY_VERSION = "1.0.0"
const _LATLON_VARIANTS = (:regular, :reduced_gaussian)
const _LATLON_POLE_POLICIES = (:none, :average, :fold)
const _LATLON_METRIC_NAMES = (
    :J,
    :g_lonlon, :g_latlat, :g_lonlat,
    :ginv_lonlon, :ginv_latlat, :ginv_lonlat,
    :area,
)

# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

struct LatLonGrid{T <: AbstractFloat} <: AbstractCurvilinearGrid
    variant::Symbol
    nlat::Int
    nlon_per_row::Vector{Int}
    R::T
    dtype::String
    ghosts::Int
    pole_policy::Symbol
    lon_start::T
    lat_edges::Vector{T}
    lat_centers::Vector{T}
end

Base.ndims(::LatLonGrid) = 2
Base.eltype(::LatLonGrid{T}) where {T} = T

n_cells(g::LatLonGrid) = sum(g.nlon_per_row)
n_vertices(g::LatLonGrid) = sum(g.nlon_per_row) + length(g.nlon_per_row)
n_edges(g::LatLonGrid) = 2 * sum(g.nlon_per_row) + length(g.nlon_per_row)
family(::LatLonGrid) = "lat_lon"

"""
    nlon(g::LatLonGrid, j::Integer) -> Int

Number of longitude cells in row `j` (1-based).
"""
function nlon(g::LatLonGrid, j::Integer)
    jj = Int(j)
    (1 ≤ jj ≤ g.nlat) ||
        throw(ArgumentError("lat_lon: j out of range [1, $(g.nlat)]: $j"))
    return g.nlon_per_row[jj]
end

"""
    nlon_uniform(g::LatLonGrid) -> Union{Int,Nothing}

Uniform longitude cell count for the regular variant; `nothing` otherwise.
"""
function nlon_uniform(g::LatLonGrid)
    g.variant === :regular || return nothing
    return g.nlon_per_row[1]
end

"""
    row_offset(g::LatLonGrid, j::Integer) -> Int

Flat starting index (0-based, row-major ragged layout) of row `j` (1-based).
"""
function row_offset(g::LatLonGrid, j::Integer)
    jj = Int(j)
    (1 ≤ jj ≤ g.nlat + 1) ||
        throw(ArgumentError("lat_lon: row offset query out of range [1, $(g.nlat + 1)]: $j"))
    off = 0
    @inbounds for k in 1:(jj - 1)
        off += g.nlon_per_row[k]
    end
    return off
end

# ---------------------------------------------------------------------------
# Parse + validation helpers
# ---------------------------------------------------------------------------

function _latlon_parse_variant(v::Symbol)
    v in _LATLON_VARIANTS ||
        throw(ArgumentError(
            "lat_lon: variant must be :regular or :reduced_gaussian; got :$v"
        ))
    return v
end
_latlon_parse_variant(v::AbstractString) = _latlon_parse_variant(Symbol(v))

function _latlon_parse_pole_policy(p::Symbol)
    p in _LATLON_POLE_POLICIES ||
        throw(ArgumentError(
            "lat_lon: pole_policy must be one of :none, :average, :fold; got :$p"
        ))
    p === :none ||
        throw(ArgumentError(
            "lat_lon: non-:none pole policies (:average, :fold) are declared but " *
                "not implemented in this phase"
        ))
    return p
end
_latlon_parse_pole_policy(p::AbstractString) = _latlon_parse_pole_policy(Symbol(p))

function _latlon_validate_lat_edges(edges::Vector{T}, nlat::Int) where {T}
    length(edges) == nlat + 1 ||
        throw(ArgumentError(
            "lat_lon: lat_edges length $(length(edges)) does not match nlat+1=$(nlat + 1)"
        ))
    for x in edges
        isfinite(x) ||
            throw(DomainError(x, "lat_lon: lat_edges must be finite"))
    end
    for k in 1:nlat
        edges[k + 1] > edges[k] ||
            throw(DomainError(edges,
                "lat_lon: lat_edges must be strictly increasing"))
    end
    tol = T(1.0e-12)
    (edges[1] ≥ -T(pi) / 2 - tol && edges[end] ≤ T(pi) / 2 + tol) ||
        throw(DomainError(edges,
            "lat_lon: lat_edges must lie within [-pi/2, pi/2]"))
    return nothing
end

function _latlon_validate_lat_centers(
        centers::Vector{T}, edges::Vector{T}, nlat::Int
    ) where {T}
    length(centers) == nlat ||
        throw(ArgumentError(
            "lat_lon: lat_centers length $(length(centers)) does not match nlat=$nlat"
        ))
    for x in centers
        isfinite(x) ||
            throw(DomainError(x, "lat_lon: lat_centers must be finite"))
    end
    for k in 1:nlat
        (edges[k] ≤ centers[k] ≤ edges[k + 1]) ||
            throw(DomainError(centers[k],
                "lat_lon: lat_centers[$k]=$(centers[k]) outside enclosing edges " *
                    "[$(edges[k]), $(edges[k + 1])]"))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Generator (entry point)
# ---------------------------------------------------------------------------

"""
    _latlon(; variant=:regular, nlon=nothing, nlat=nothing,
              nlon_per_row=nothing, lat_edges=nothing, lat_centers=nothing,
              R=6.371e6, dtype=Float64, ghosts=0, pole_policy=:none,
              lon_start=nothing) -> LatLonGrid{T}

Primary constructor for the lat-lon family. See `docs/GRIDS_API.md` §2.3 for
the cross-binding signature contract. Public entry point (Julia):
`EarthSciDiscretizations.grids.lat_lon`.
"""
function _latlon(;
        variant = :regular,
        nlon::Union{Int, Nothing} = nothing,
        nlat::Union{Int, Nothing} = nothing,
        nlon_per_row = nothing,
        lat_edges = nothing,
        lat_centers = nothing,
        R::Real = 6.371e6,
        dtype::Type{T} = Float64,
        ghosts::Int = 0,
        pole_policy = :none,
        lon_start = nothing
    ) where {T <: AbstractFloat}

    var = _latlon_parse_variant(variant isa Symbol ? variant : Symbol(variant))
    pol = _latlon_parse_pole_policy(pole_policy isa Symbol ? pole_policy : Symbol(pole_policy))

    ghosts ≥ 0 ||
        throw(ArgumentError("lat_lon: ghosts must be ≥ 0; got $ghosts"))

    R_T = T(R)
    (isfinite(R_T) && R_T > zero(T)) ||
        throw(DomainError(R, "lat_lon: R must be positive and finite"))

    lon_start_T = lon_start === nothing ? -T(pi) : T(lon_start)
    isfinite(lon_start_T) ||
        throw(DomainError(lon_start, "lat_lon: lon_start must be finite"))

    local nlat_i::Int
    local per_row::Vector{Int}

    if var === :regular
        nlon_per_row === nothing ||
            throw(ArgumentError("lat_lon: nlon_per_row is not allowed for variant=:regular"))
        nlon === nothing &&
            throw(ArgumentError("lat_lon: missing required keyword argument `nlon`"))
        nlat === nothing &&
            throw(ArgumentError("lat_lon: missing required keyword argument `nlat`"))
        nlon ≥ 1 ||
            throw(DomainError(nlon, "lat_lon: nlon must be ≥ 1"))
        nlat ≥ 1 ||
            throw(DomainError(nlat, "lat_lon: nlat must be ≥ 1"))
        nlat_i = nlat
        per_row = fill(Int(nlon), nlat_i)
    else
        # :reduced_gaussian
        nlon === nothing ||
            throw(ArgumentError(
                "lat_lon: nlon is not allowed for variant=:reduced_gaussian; " *
                    "use nlon_per_row"
            ))
        nlon_per_row === nothing &&
            throw(ArgumentError(
                "lat_lon: missing required keyword argument `nlon_per_row` " *
                    "for variant=:reduced_gaussian"
            ))
        per_row = Int[Int(x) for x in nlon_per_row]
        nlat_i = nlat === nothing ? length(per_row) : Int(nlat)
        nlat_i ≥ 1 ||
            throw(DomainError(nlat_i, "lat_lon: nlat must be ≥ 1"))
        length(per_row) == nlat_i ||
            throw(ArgumentError(
                "lat_lon: nlon_per_row length $(length(per_row)) does not match " *
                    "nlat=$nlat_i"
            ))
        for (j, n) in enumerate(per_row)
            n ≥ 1 ||
                throw(DomainError(n, "lat_lon: nlon_per_row[$j]=$n must be ≥ 1"))
        end
    end

    if lat_edges === nothing
        dlat = T(pi) / T(nlat_i)
        edges = T[-T(pi) / 2 + T(k) * dlat for k in 0:nlat_i]
    else
        edges = T[T(x) for x in lat_edges]
        _latlon_validate_lat_edges(edges, nlat_i)
    end

    centers = if lat_centers === nothing
        T[(edges[k] + edges[k + 1]) / T(2) for k in 1:nlat_i]
    else
        c = T[T(x) for x in lat_centers]
        _latlon_validate_lat_centers(c, edges, nlat_i)
        c
    end

    return LatLonGrid{T}(
        var, nlat_i, per_row, R_T, _dtype_string(T), ghosts,
        pol, lon_start_T, edges, centers,
    )
end

# ---------------------------------------------------------------------------
# Bounds helpers
# ---------------------------------------------------------------------------

function _latlon_check_cell(g::LatLonGrid, j::Integer, i::Integer)
    jj = Int(j); ii = Int(i)
    (1 ≤ jj ≤ g.nlat) ||
        throw(ArgumentError("lat_lon: j out of range [1, $(g.nlat)]: $j"))
    n = g.nlon_per_row[jj]
    (1 ≤ ii ≤ n) ||
        throw(ArgumentError("lat_lon: i out of range [1, $n] for row $jj: $i"))
    return nothing
end

# Map 1-based column ii in a row of width from_n to the nearest-center 1-based
# column in a row of width to_n. Mirrors the Python helper `_map_i` modulo the
# 0-/1-based index convention.
function _latlon_map_i(ii::Int, from_n::Int, to_n::Int)
    from_n == to_n && return ii
    frac = (Float64(ii) - 0.5) / Float64(from_n)
    k = Int(floor(frac * Float64(to_n)))
    return clamp(k, 0, to_n - 1) + 1
end

# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------

"""
    cell_centers(g::LatLonGrid, j::Integer, i::Integer) -> (lon, lat)

Cell-centre geographic coordinates (radians) for cell `(j, i)` (1-based:
`1 ≤ j ≤ nlat`, `1 ≤ i ≤ nlon_per_row[j]`).
"""
function cell_centers(g::LatLonGrid{T}, j::Integer, i::Integer) where {T}
    _latlon_check_cell(g, j, i)
    jj = Int(j); ii = Int(i)
    n = g.nlon_per_row[jj]
    dlon = T(2) * T(pi) / T(n)
    lon = g.lon_start + (T(ii) - T(0.5)) * dlon
    return (lon, g.lat_centers[jj])
end

"""
    cell_centers(g::LatLonGrid) -> NamedTuple{(:lon,:lat)}

Flat ragged-row-major `(lon, lat)` arrays of length `n_cells(g)`. Row `j`
occupies indices `row_offset(g, j) + 1 : row_offset(g, j) + nlon(g, j)`.
"""
function cell_centers(g::LatLonGrid{T}) where {T}
    nc = n_cells(g)
    lon = Vector{T}(undef, nc)
    lat = Vector{T}(undef, nc)
    k = 0
    @inbounds for j in 1:g.nlat
        n = g.nlon_per_row[j]
        dlon = T(2) * T(pi) / T(n)
        lat_c = g.lat_centers[j]
        for i in 1:n
            k += 1
            lon[k] = g.lon_start + (T(i) - T(0.5)) * dlon
            lat[k] = lat_c
        end
    end
    return (lon = lon, lat = lat)
end

"""
    lon_edges(g::LatLonGrid, j::Integer) -> Vector{T}

Longitude cell edges for row `j` (length `nlon_per_row[j] + 1`).
"""
function lon_edges(g::LatLonGrid{T}, j::Integer) where {T}
    jj = Int(j)
    (1 ≤ jj ≤ g.nlat) ||
        throw(ArgumentError("lat_lon: j out of range [1, $(g.nlat)]: $j"))
    n = g.nlon_per_row[jj]
    dlon = T(2) * T(pi) / T(n)
    return T[g.lon_start + T(k) * dlon for k in 0:n]
end

"""
    lon_centers(g::LatLonGrid, j::Integer) -> Vector{T}

Longitude cell centres for row `j` (length `nlon_per_row[j]`).
"""
function lon_centers(g::LatLonGrid{T}, j::Integer) where {T}
    jj = Int(j)
    (1 ≤ jj ≤ g.nlat) ||
        throw(ArgumentError("lat_lon: j out of range [1, $(g.nlat)]: $j"))
    n = g.nlon_per_row[jj]
    dlon = T(2) * T(pi) / T(n)
    return T[g.lon_start + (T(k) - T(0.5)) * dlon for k in 1:n]
end

"""
    neighbors(g::LatLonGrid, j::Integer, i::Integer)
        -> Dict{Symbol, Union{Tuple{Int,Int}, Nothing}}

Face neighbours of cell `(j, i)` keyed by `:W`, `:E`, `:S`, `:N`. Longitude
wraps periodically. Under the default `pole_policy=:none`, the `:S` neighbour
of the first row and the `:N` neighbour of the last row are `nothing`. For
reduced-Gaussian grids the N/S neighbour is the nearest-centre cell in the
adjacent row (accounting for differing `nlon`).
"""
function neighbors(g::LatLonGrid, j::Integer, i::Integer)
    _latlon_check_cell(g, j, i)
    jj = Int(j); ii = Int(i)
    n_i = g.nlon_per_row[jj]
    out = Dict{Symbol, Union{Tuple{Int, Int}, Nothing}}()
    out[:W] = (jj, ii == 1 ? n_i : ii - 1)
    out[:E] = (jj, ii == n_i ? 1 : ii + 1)
    if jj == 1
        out[:S] = _latlon_pole_neighbor(g)
    else
        n_s = g.nlon_per_row[jj - 1]
        out[:S] = (jj - 1, _latlon_map_i(ii, n_i, n_s))
    end
    if jj == g.nlat
        out[:N] = _latlon_pole_neighbor(g)
    else
        n_n = g.nlon_per_row[jj + 1]
        out[:N] = (jj + 1, _latlon_map_i(ii, n_i, n_n))
    end
    return out
end

_latlon_pole_neighbor(::LatLonGrid) = nothing  # :none is the only implemented policy

"""
    cell_area(g::LatLonGrid, j::Integer, i::Integer) -> T

Spherical-rectangle area for cell `(j, i)`:
`R² · dlon · (sin(lat_n) − sin(lat_s))`.
"""
function cell_area(g::LatLonGrid{T}, j::Integer, i::Integer) where {T}
    _latlon_check_cell(g, j, i)
    jj = Int(j)
    n = g.nlon_per_row[jj]
    dlon = T(2) * T(pi) / T(n)
    lat_s = g.lat_edges[jj]
    lat_n = g.lat_edges[jj + 1]
    return g.R * g.R * dlon * (sin(lat_n) - sin(lat_s))
end

"""
    metric_eval(g::LatLonGrid, name::Symbol, j::Integer, i::Integer) -> T

Evaluate a named metric at cell centre `(j, i)`. Supported names:

- `:J`            — Jacobian determinant `R² |cos(lat)|`
- `:g_lonlon`, `:g_latlat`, `:g_lonlat`  — covariant metric tensor components
- `:ginv_lonlon`, `:ginv_latlat`, `:ginv_lonlat` — inverse metric
- `:area`         — spherical-rectangle cell area

The lat-lon metric is longitudinally independent, so `i` is unused for the
non-area metrics (bounds are still checked).
"""
function metric_eval(g::LatLonGrid{T}, name::Symbol, j::Integer, i::Integer) where {T}
    name in _LATLON_METRIC_NAMES ||
        throw(ArgumentError("lat_lon: unknown metric :$name"))
    _latlon_check_cell(g, j, i)
    name === :area && return cell_area(g, j, i)
    lat = g.lat_centers[Int(j)]
    cos_lat = cos(lat)
    r2 = g.R * g.R
    g_ll = r2 * cos_lat * cos_lat
    g_pp = r2
    if name === :J
        return r2 * abs(cos_lat)
    elseif name === :g_lonlon
        return g_ll
    elseif name === :g_latlat
        return g_pp
    elseif name === :g_lonlat
        return zero(T)
    elseif name === :ginv_lonlon
        return g_ll > zero(T) ? one(T) / g_ll : T(Inf)
    elseif name === :ginv_latlat
        return one(T) / g_pp
    else  # :ginv_lonlat
        return zero(T)
    end
end

metric_eval(g::LatLonGrid, name::AbstractString, j::Integer, i::Integer) =
    metric_eval(g, Symbol(name), j, i)

# ---------------------------------------------------------------------------
# Provenance + .esm lowering
# ---------------------------------------------------------------------------

function _latlon_generator_name(g::LatLonGrid)
    return g.variant === :regular ? "lat_lon_regular" : "lat_lon_reduced_gaussian"
end

function _latlon_provenance(g::LatLonGrid)
    return Dict{String, Any}(
        "binding" => "julia",
        "binding_version" => _LATLON_FAMILY_VERSION,
        "source" => "EarthSciDiscretizations.grids.lat_lon",
        "generator" => _latlon_generator_name(g),
    )
end

"""
    to_esm(g::LatLonGrid) -> Dict{String,Any}

Declarative `.esm` lowering per the 2026-04-20 scope correction. The output
carries `family`, `version`, `dtype`, `topology`, `variant`, `generator`,
`params`, and `provenance`; derived per-cell geometry never appears in the wire
form. Matches the Python reference binding's key set so cross-binding byte
equality is achievable by the conformance harness.
"""
function to_esm(g::LatLonGrid)
    params = if g.variant === :regular
        Dict{String, Any}(
            "nlon" => Int(g.nlon_per_row[1]),
            "nlat" => Int(g.nlat),
            "R" => Float64(g.R),
            "ghosts" => Int(g.ghosts),
            "pole_policy" => String(g.pole_policy),
            "lon_start" => Float64(g.lon_start),
        )
    else
        Dict{String, Any}(
            "nlat" => Int(g.nlat),
            "nlon_per_row" => Int[x for x in g.nlon_per_row],
            "lat_edges" => Float64[Float64(x) for x in g.lat_edges],
            "R" => Float64(g.R),
            "ghosts" => Int(g.ghosts),
            "pole_policy" => String(g.pole_policy),
            "lon_start" => Float64(g.lon_start),
        )
    end
    return Dict{String, Any}(
        "family" => "lat_lon",
        "version" => _LATLON_FAMILY_VERSION,
        "dtype" => g.dtype,
        "topology" => "rectilinear",
        "variant" => String(g.variant),
        "generator" => _latlon_generator_name(g),
        "params" => params,
        "provenance" => _latlon_provenance(g),
    )
end

# Public namespace registration: the `lat_lon` entry point in the
# `EarthSciDiscretizations.grids` submodule (defined at the bottom of
# `src/EarthSciDiscretizations.jl`) aliases `_latlon` here. See
# `GRIDS_API.md` §2.3 for the public call form.

# ---------------------------------------------------------------------------
# ESS Grid trait — Tier C + Tier M (sphere-surface curvilinear)
#
# Flat ragged-row-major cell layout: cell `(j, i)` (1-based) lives at flat
# index `row_offset(g, j) + i`. Longitude wraps periodically. Latitude
# stepping uses `_latlon_map_i` to handle reduced-Gaussian rows where
# `nlon_per_row[j]` differs across `j`. Pole neighbours are `0` under the
# default `pole_policy=:none`.
# ---------------------------------------------------------------------------

n_dims(::LatLonGrid) = 2
axis_names(::LatLonGrid) = (:lon, :lat)

function _latlon_axis_check(::LatLonGrid, axis::Symbol)
    axis in (:lon, :lat) ||
        throw(ArgumentError("lat_lon: unknown axis :$axis (expected :lon or :lat)"))
    return nothing
end

function cell_centers(g::LatLonGrid{T}, axis::Symbol) where {T}
    _latlon_axis_check(g, axis)
    return _grid_memo!(g, (:cell_centers, axis)) do
        nt = cell_centers(g)  # NamedTuple{(:lon,:lat)}
        return axis === :lon ? nt.lon : nt.lat
    end
end

function cell_widths(g::LatLonGrid{T}, axis::Symbol) where {T}
    _latlon_axis_check(g, axis)
    return _grid_memo!(g, (:cell_widths, axis)) do
        nc = n_cells(g)
        out = Vector{T}(undef, nc)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            if axis === :lon
                w = T(2) * T(pi) / T(n_i)
                for _ in 1:n_i
                    k += 1
                    out[k] = w
                end
            else
                w = g.lat_edges[j + 1] - g.lat_edges[j]
                for _ in 1:n_i
                    k += 1
                    out[k] = w
                end
            end
        end
        return out
    end
end

function cell_volume(g::LatLonGrid{T}) where {T}
    return _grid_memo!(g, :cell_volume) do
        nc = n_cells(g)
        out = Vector{T}(undef, nc)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            dlon = T(2) * T(pi) / T(n_i)
            sin_n = sin(g.lat_edges[j + 1])
            sin_s = sin(g.lat_edges[j])
            area = g.R * g.R * dlon * (sin_n - sin_s)
            for _ in 1:n_i
                k += 1
                out[k] = area
            end
        end
        return out
    end
end

function neighbor_indices(g::LatLonGrid, axis::Symbol, offset::Int)
    _latlon_axis_check(g, axis)
    return _grid_memo!(g, (:neighbor_indices, axis, offset)) do
        nc = n_cells(g)
        out = Vector{Int}(undef, nc)
        # Precompute row offsets once.
        row_off = Vector{Int}(undef, g.nlat + 1)
        row_off[1] = 0
        @inbounds for j in 1:g.nlat
            row_off[j + 1] = row_off[j] + g.nlon_per_row[j]
        end
        if axis === :lon
            @inbounds for j in 1:g.nlat
                n_i = g.nlon_per_row[j]
                base = row_off[j]
                for i in 1:n_i
                    ii = mod(i - 1 + offset, n_i) + 1
                    out[base + i] = row_off[j] + ii
                end
            end
        else  # :lat
            @inbounds for j in 1:g.nlat
                n_i = g.nlon_per_row[j]
                base = row_off[j]
                jj = j + offset
                if 1 <= jj <= g.nlat
                    n_t = g.nlon_per_row[jj]
                    base_t = row_off[jj]
                    for i in 1:n_i
                        i_t = _latlon_map_i(i, n_i, n_t)
                        out[base + i] = base_t + i_t
                    end
                else
                    # Pole: pole_policy=:none → no neighbour.
                    for i in 1:n_i
                        out[base + i] = 0
                    end
                end
            end
        end
        return out
    end
end

function boundary_mask(g::LatLonGrid, axis::Symbol, side::Symbol)
    _latlon_axis_check(g, axis)
    side in (:lower, :upper) ||
        throw(ArgumentError("lat_lon: side must be :lower or :upper; got :$side"))
    return _grid_memo!(g, (:boundary_mask, axis, side)) do
        nc = n_cells(g)
        out = falses(nc)
        if axis === :lon
            # Longitude wraps periodically — no lon boundary cells.
            return out
        end
        target = side === :lower ? 1 : g.nlat
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            mark = j == target
            for _ in 1:n_i
                k += 1
                out[k] = mark
            end
        end
        return out
    end
end

# Tier-M curvilinear metric. The lat-lon basis on a sphere of radius R has
# `g_λλ = R² cos²(φ)`, `g_φφ = R²`, `g_λφ = 0`. The longitudinally-uniform
# layout means every cell in row `j` shares the same metric tensor.

function metric_g(g::LatLonGrid{T}) where {T}
    return _grid_memo!(g, :metric_g) do
        nc = n_cells(g)
        out = zeros(T, nc, 2, 2)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            cos_lat = cos(g.lat_centers[j])
            g_ll = g.R * g.R * cos_lat * cos_lat
            g_pp = g.R * g.R
            for _ in 1:n_i
                k += 1
                out[k, 1, 1] = g_ll
                out[k, 2, 2] = g_pp
            end
        end
        return out
    end
end

function metric_ginv(g::LatLonGrid{T}) where {T}
    return _grid_memo!(g, :metric_ginv) do
        nc = n_cells(g)
        out = zeros(T, nc, 2, 2)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            cos_lat = cos(g.lat_centers[j])
            g_ll = g.R * g.R * cos_lat * cos_lat
            g_pp = g.R * g.R
            inv_ll = g_ll > zero(T) ? one(T) / g_ll : T(Inf)
            inv_pp = one(T) / g_pp
            for _ in 1:n_i
                k += 1
                out[k, 1, 1] = inv_ll
                out[k, 2, 2] = inv_pp
            end
        end
        return out
    end
end

function metric_jacobian(g::LatLonGrid{T}) where {T}
    return _grid_memo!(g, :metric_jacobian) do
        nc = n_cells(g)
        out = Vector{T}(undef, nc)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            J = g.R * g.R * abs(cos(g.lat_centers[j]))
            for _ in 1:n_i
                k += 1
                out[k] = J
            end
        end
        return out
    end
end

function metric_dgij_dxk(g::LatLonGrid{T}) where {T}
    # Lat-lon layout uses physical (lon, lat) as the computational axes too,
    # so `∂g_ij/∂x^k` reduces to derivatives w.r.t. (lon, lat). All entries
    # vanish except `∂g_λλ/∂φ = -2 R² cos(φ) sin(φ)`.
    return _grid_memo!(g, :metric_dgij_dxk) do
        nc = n_cells(g)
        out = zeros(T, nc, 2, 2, 2)
        k = 0
        @inbounds for j in 1:g.nlat
            n_i = g.nlon_per_row[j]
            φ = g.lat_centers[j]
            dgll_dφ = -T(2) * g.R * g.R * cos(φ) * sin(φ)
            for _ in 1:n_i
                k += 1
                out[k, 1, 1, 2] = dgll_dφ  # ∂g_λλ / ∂φ
            end
        end
        return out
    end
end

function coord_jacobian(g::LatLonGrid{T}, target::Symbol) where {T}
    target === :lon_lat ||
        throw(ArgumentError("lat_lon: coord_jacobian only supports target=:lon_lat; got :$target"))
    # Computational and target axes coincide, so `∂(comp)/∂(target) = δ`.
    return _grid_memo!(g, (:coord_jacobian, target)) do
        nc = n_cells(g)
        out = zeros(T, nc, 2, 2)
        @inbounds for k in 1:nc
            out[k, 1, 1] = one(T)
            out[k, 2, 2] = one(T)
        end
        return out
    end
end

function coord_jacobian_second(g::LatLonGrid{T}, target::Symbol) where {T}
    target === :lon_lat ||
        throw(ArgumentError("lat_lon: coord_jacobian_second only supports target=:lon_lat; got :$target"))
    return _grid_memo!(g, (:coord_jacobian_second, target)) do
        zeros(T, n_cells(g), 2, 2, 2)
    end
end
