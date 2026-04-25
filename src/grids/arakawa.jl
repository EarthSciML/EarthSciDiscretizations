"""
Arakawa staggering accessor runtime.

Arakawa staggering is a *transform* over an underlying base grid. Given a base
grid (Cartesian, lat-lon, …) and a stagger label (A/B/C/D/E), the arakawa
runtime provides on-demand accessors for cell-center, u-face, v-face, and
corner locations per the staggering convention.

Per GRIDS_API.md §0 correction: the `.esm` lowering is a small declarative
config (family + base-grid ref + stagger + dimensions + extents). Geometry
is derived from that config by the accessors here, not serialized as a blob.

Staggering conventions (2-D, horizontal):

- A: h, u, v colocated at cell centers.
- B: h at cell centers; u, v colocated at corners.
- C: h at cell centers; u at u-faces (east-west); v at v-faces (north-south).
- D: h at cell centers; u at v-faces; v at u-faces (swapped from C).
- E: rotated B-grid; topologically equivalent to B with a 45° rotation flag
  in the lowered config.
"""

@enum ArakawaStagger::UInt8 ArakawaA ArakawaB ArakawaC ArakawaD ArakawaE

# Reuse VarLocation from staggering.jl (loaded before this file via
# EarthSciDiscretizations.jl include order).

# Variable-to-location table per stagger. Returns (h_loc, u_loc, v_loc).
function arakawa_variable_locations(s::ArakawaStagger)
    return s === ArakawaA ? (CellCenter, CellCenter, CellCenter) :
        s === ArakawaB ? (CellCenter, Corner, Corner) :
        s === ArakawaC ? (CellCenter, UEdge, VEdge) :
        s === ArakawaD ? (CellCenter, VEdge, UEdge) :
        (CellCenter, Corner, Corner)  # ArakawaE: like B, rotated 45°
end

# Location shape in (ni, nj) on a base grid of (nx, ny) interior cells.
function arakawa_location_shape(loc::VarLocation, nx::Int, ny::Int)
    return loc === CellCenter ? (nx, ny) :
        loc === UEdge ? (nx + 1, ny) :
        loc === VEdge ? (nx, ny + 1) :
        (nx + 1, ny + 1)   # Corner
end

# ---------------------------------------------------------------------------
# Base grid abstraction
#
# Phase 1 (cartesian) and Phase 2/3 (lat-lon) grid types are not yet on main.
# The arakawa runtime is parameterised over an `ArakawaBaseGrid`; a minimal
# `CartesianBase` is provided here so the accessors are exercisable today.
# Future Phase 1/2/3 grid structs can subtype `ArakawaBaseGrid` and supply
# the same four primitives (cell-center, x-edge, y-edge, corner coordinates
# plus cell counts) to plug in without touching this file.
# ---------------------------------------------------------------------------

abstract type ArakawaBaseGrid end

struct CartesianBase{T} <: ArakawaBaseGrid
    xlo::T
    xhi::T
    ylo::T
    yhi::T
    nx::Int
    ny::Int
end

function CartesianBase(; xlo::Real, xhi::Real, ylo::Real, yhi::Real, nx::Int, ny::Int)
    nx > 0 || throw(DomainError(nx, "nx must be positive"))
    ny > 0 || throw(DomainError(ny, "ny must be positive"))
    xhi > xlo || throw(DomainError(xhi, "xhi must be greater than xlo"))
    yhi > ylo || throw(DomainError(yhi, "yhi must be greater than ylo"))
    T = promote_type(typeof(xlo), typeof(xhi), typeof(ylo), typeof(yhi))
    return CartesianBase{T}(T(xlo), T(xhi), T(ylo), T(yhi), nx, ny)
end

arakawa_nx(b::CartesianBase) = b.nx
arakawa_ny(b::CartesianBase) = b.ny
arakawa_dx(b::CartesianBase{T}) where {T} = (b.xhi - b.xlo) / T(b.nx)
arakawa_dy(b::CartesianBase{T}) where {T} = (b.yhi - b.ylo) / T(b.ny)

function arakawa_cell_center(b::CartesianBase{T}, i::Int, j::Int) where {T}
    return (
        b.xlo + (T(i) - T(0.5)) * arakawa_dx(b),
        b.ylo + (T(j) - T(0.5)) * arakawa_dy(b),
    )
end

function arakawa_x_edge(b::CartesianBase{T}, i::Int, j::Int) where {T}
    return (
        b.xlo + T(i - 1) * arakawa_dx(b),
        b.ylo + (T(j) - T(0.5)) * arakawa_dy(b),
    )
end

function arakawa_y_edge(b::CartesianBase{T}, i::Int, j::Int) where {T}
    return (
        b.xlo + (T(i) - T(0.5)) * arakawa_dx(b),
        b.ylo + T(j - 1) * arakawa_dy(b),
    )
end

function arakawa_corner(b::CartesianBase{T}, i::Int, j::Int) where {T}
    return (
        b.xlo + T(i - 1) * arakawa_dx(b),
        b.ylo + T(j - 1) * arakawa_dy(b),
    )
end

_arakawa_base_esm(b::CartesianBase) = Dict{String, Any}(
    "family" => "cartesian",
    "nx" => b.nx,
    "ny" => b.ny,
    "extent" => [[b.xlo, b.ylo], [b.xhi, b.yhi]],
)

# ---------------------------------------------------------------------------
# ArakawaGrid: a base grid + stagger label
# ---------------------------------------------------------------------------

struct ArakawaGrid{T, B <: ArakawaBaseGrid} <: AbstractStaggeredGrid
    base::B
    stagger::ArakawaStagger
    ghosts::Int
end

function ArakawaGrid(
        base::B, stagger::ArakawaStagger;
        ghosts::Int = 0, dtype::Type{T} = Float64
    ) where {T, B <: ArakawaBaseGrid}
    ghosts >= 0 || throw(DomainError(ghosts, "ghosts must be non-negative"))
    return ArakawaGrid{T, B}(base, stagger, ghosts)
end

Base.eltype(::Type{ArakawaGrid{T, B}}) where {T, B} = T
Base.eltype(g::ArakawaGrid) = eltype(typeof(g))

arakawa_nx(g::ArakawaGrid) = arakawa_nx(g.base)
arakawa_ny(g::ArakawaGrid) = arakawa_ny(g.base)

# Shape of a specific location on this grid (interior cells only; ghosts
# are reported separately via `ghosts` and are not expanded here).
arakawa_shape(g::ArakawaGrid, loc::VarLocation) =
    arakawa_location_shape(loc, arakawa_nx(g), arakawa_ny(g))

# Shape of the location where each named variable lives (h, u, v).
function variable_shape(g::ArakawaGrid, var::Symbol)
    h_loc, u_loc, v_loc = arakawa_variable_locations(g.stagger)
    loc = var === :h ? h_loc :
        var === :u ? u_loc :
        var === :v ? v_loc :
        throw(ArgumentError("Unknown arakawa variable $var; expected :h, :u, or :v"))
    return arakawa_shape(g, loc)
end

# ---------------------------------------------------------------------------
# Accessors — DERIVE geometry from the declaration (mayor's correction).
# ---------------------------------------------------------------------------

"""
    cell_centers(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of cell-centre `(i, j)` (`1`-based). `(i, j)` must satisfy
`1 ≤ i ≤ nx` and `1 ≤ j ≤ ny`.
"""
function cell_centers(g::ArakawaGrid, i::Int, j::Int)
    _check_bounds(g, CellCenter, i, j)
    return arakawa_cell_center(g.base, i, j)
end

"""
    u_face(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of the location where `u` lives at index `(i, j)` for this
stagger. Valid index range depends on the stagger (see `variable_shape`).
"""
function u_face(g::ArakawaGrid, i::Int, j::Int)
    _, u_loc, _ = arakawa_variable_locations(g.stagger)
    _check_bounds(g, u_loc, i, j)
    return _coord_at(g, u_loc, i, j)
end

"""
    v_face(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of the location where `v` lives at index `(i, j)` for this
stagger.
"""
function v_face(g::ArakawaGrid, i::Int, j::Int)
    _, _, v_loc = arakawa_variable_locations(g.stagger)
    _check_bounds(g, v_loc, i, j)
    return _coord_at(g, v_loc, i, j)
end

"""
    corners(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of cell-corner `(i, j)` (`1 ≤ i ≤ nx+1`, `1 ≤ j ≤ ny+1`).
"""
function corners(g::ArakawaGrid, i::Int, j::Int)
    _check_bounds(g, Corner, i, j)
    return arakawa_corner(g.base, i, j)
end

function _coord_at(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    return loc === CellCenter ? arakawa_cell_center(g.base, i, j) :
        loc === UEdge ? arakawa_x_edge(g.base, i, j) :
        loc === VEdge ? arakawa_y_edge(g.base, i, j) :
        arakawa_corner(g.base, i, j)
end

function _check_bounds(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    ni, nj = arakawa_shape(g, loc)
    return (1 <= i <= ni && 1 <= j <= nj) ||
        throw(DomainError((i, j), "index out of bounds for location $loc with shape ($ni, $nj)"))
end

"""
    neighbors(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
        -> (west, east, south, north)

Return the four axial neighbours of `(loc, i, j)` as `(loc, i, j)` tuples, or
`nothing` at domain boundaries. Does not cross ghost cells (ghost halos are
reported via `g.ghosts` but not expanded into the index range here).
"""
function neighbors(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    _check_bounds(g, loc, i, j)
    ni, nj = arakawa_shape(g, loc)
    west = i > 1 ? (loc, i - 1, j) : nothing
    east = i < ni ? (loc, i + 1, j) : nothing
    south = j > 1 ? (loc, i, j - 1) : nothing
    north = j < nj ? (loc, i, j + 1) : nothing
    return (west, east, south, north)
end

neighbors(g::ArakawaGrid, i::Int, j::Int) = neighbors(g, CellCenter, i, j)

"""
    metric_eval(g::ArakawaGrid, name::Symbol, i::Int, j::Int) -> scalar

Return the named metric at cell `(i, j)`. Supported names:

- `:dx`   — zonal spacing at `(i, j)`
- `:dy`   — meridional spacing at `(i, j)`
- `:area` — cell-centre area / volume

For `CartesianBase`, metrics are uniform and the `(i, j)` arguments are
ignored (bounds are still checked).
"""
function metric_eval(g::ArakawaGrid{T}, name::Symbol, i::Int, j::Int) where {T}
    _check_bounds(g, CellCenter, i, j)
    dx = T(arakawa_dx(g.base))
    dy = T(arakawa_dy(g.base))
    return if name === :dx
        dx
    elseif name === :dy
        dy
    elseif name === :area
        dx * dy
    else
        throw(ArgumentError("Unknown metric $name; expected :dx, :dy, or :area"))
    end
end

# ---------------------------------------------------------------------------
# .esm lowering: declarative config per §6 schema + the §0 correction.
# ---------------------------------------------------------------------------

"""
    to_esm(g::ArakawaGrid) -> Dict{String,Any}

Lower the grid to a §6-schema-valid `.esm` declaration. The result is a
small config — `family`, `stagger`, `base`, `dtype`, `ghosts`, `n_cells`,
`topology` — without any inline geometry arrays.
"""
function to_esm(g::ArakawaGrid{T}) where {T}
    return Dict{String, Any}(
        "family" => "arakawa",
        "dtype" => _dtype_string(T),
        "topology" => "block_structured",
        "ghosts" => g.ghosts,
        "n_cells" => arakawa_nx(g) * arakawa_ny(g),
        "stagger" => _stagger_string(g.stagger),
        "rotated" => g.stagger === ArakawaE,
        "base" => _arakawa_base_esm(g.base),
    )
end

function _stagger_string(s::ArakawaStagger)
    return s === ArakawaA ? "A" :
        s === ArakawaB ? "B" :
        s === ArakawaC ? "C" :
        s === ArakawaD ? "D" : "E"
end

_dtype_string(::Type{Float64}) = "float64"
_dtype_string(::Type{Float32}) = "float32"
_dtype_string(::Type{T}) where {T} = string(T)

# ---------------------------------------------------------------------------
# ESS Grid trait — Tier C + Tier S (staggered horizontal mesh)
#
# Cell-centered Tier-C arrays are flat row-major over the (nx, ny) interior
# (slow axis = j). Face / corner accessors stay scalar — they are
# stagger-specific and live above as `u_face`, `v_face`, `corners`.
# ---------------------------------------------------------------------------

n_dims(::ArakawaGrid) = 2
axis_names(::ArakawaGrid) = (:x, :y)
n_cells(g::ArakawaGrid) = arakawa_nx(g) * arakawa_ny(g)

function _arakawa_axis_idx(::ArakawaGrid, axis::Symbol)
    axis === :x && return 1
    axis === :y && return 2
    throw(ArgumentError("arakawa: unknown axis :$axis (expected :x or :y)"))
end

function cell_centers(g::ArakawaGrid{T}, axis::Symbol) where {T}
    d = _arakawa_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_centers, axis)) do
        nx = arakawa_nx(g); ny = arakawa_ny(g)
        out = Vector{T}(undef, nx * ny)
        @inbounds for j in 1:ny, i in 1:nx
            xy = arakawa_cell_center(g.base, i, j)
            out[(j - 1) * nx + i] = T(xy[d])
        end
        return out
    end
end

function cell_widths(g::ArakawaGrid{T}, axis::Symbol) where {T}
    d = _arakawa_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_widths, axis)) do
        nx = arakawa_nx(g); ny = arakawa_ny(g)
        w = d == 1 ? T(arakawa_dx(g.base)) : T(arakawa_dy(g.base))
        return fill(w, nx * ny)
    end
end

function cell_volume(g::ArakawaGrid{T}) where {T}
    return _grid_memo!(g, :cell_volume) do
        nx = arakawa_nx(g); ny = arakawa_ny(g)
        a = T(arakawa_dx(g.base)) * T(arakawa_dy(g.base))
        return fill(a, nx * ny)
    end
end

function neighbor_indices(g::ArakawaGrid, axis::Symbol, offset::Int)
    d = _arakawa_axis_idx(g, axis)
    return _grid_memo!(g, (:neighbor_indices, axis, offset)) do
        nx = arakawa_nx(g); ny = arakawa_ny(g)
        out = Vector{Int}(undef, nx * ny)
        @inbounds for j in 1:ny, i in 1:nx
            ii = d == 1 ? i + offset : i
            jj = d == 2 ? j + offset : j
            k = (j - 1) * nx + i
            out[k] = (1 <= ii <= nx && 1 <= jj <= ny) ?
                ((jj - 1) * nx + ii) : 0
        end
        return out
    end
end

function boundary_mask(g::ArakawaGrid, axis::Symbol, side::Symbol)
    d = _arakawa_axis_idx(g, axis)
    side in (:lower, :upper) ||
        throw(ArgumentError("arakawa: side must be :lower or :upper; got :$side"))
    return _grid_memo!(g, (:boundary_mask, axis, side)) do
        nx = arakawa_nx(g); ny = arakawa_ny(g)
        out = falses(nx * ny)
        target = side === :lower ? 1 : (d == 1 ? nx : ny)
        @inbounds for j in 1:ny, i in 1:nx
            v = d == 1 ? i : j
            out[(j - 1) * nx + i] = v == target
        end
        return out
    end
end
