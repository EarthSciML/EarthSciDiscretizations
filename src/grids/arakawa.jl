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
    s === ArakawaA ? (CellCenter, CellCenter, CellCenter) :
    s === ArakawaB ? (CellCenter, Corner, Corner) :
    s === ArakawaC ? (CellCenter, UEdge, VEdge) :
    s === ArakawaD ? (CellCenter, VEdge, UEdge) :
    (CellCenter, Corner, Corner)  # ArakawaE: like B, rotated 45°
end

# Location shape in (ni, nj) on a base grid of (nx, ny) interior cells.
function arakawa_location_shape(loc::VarLocation, nx::Int, ny::Int)
    loc === CellCenter ? (nx, ny) :
    loc === UEdge      ? (nx + 1, ny) :
    loc === VEdge      ? (nx, ny + 1) :
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
    CartesianBase{T}(T(xlo), T(xhi), T(ylo), T(yhi), nx, ny)
end

arakawa_nx(b::CartesianBase) = b.nx
arakawa_ny(b::CartesianBase) = b.ny
arakawa_dx(b::CartesianBase{T}) where {T} = (b.xhi - b.xlo) / T(b.nx)
arakawa_dy(b::CartesianBase{T}) where {T} = (b.yhi - b.ylo) / T(b.ny)

function arakawa_cell_center(b::CartesianBase{T}, i::Int, j::Int) where {T}
    (b.xlo + (T(i) - T(0.5)) * arakawa_dx(b),
     b.ylo + (T(j) - T(0.5)) * arakawa_dy(b))
end

function arakawa_x_edge(b::CartesianBase{T}, i::Int, j::Int) where {T}
    (b.xlo + T(i - 1) * arakawa_dx(b),
     b.ylo + (T(j) - T(0.5)) * arakawa_dy(b))
end

function arakawa_y_edge(b::CartesianBase{T}, i::Int, j::Int) where {T}
    (b.xlo + (T(i) - T(0.5)) * arakawa_dx(b),
     b.ylo + T(j - 1) * arakawa_dy(b))
end

function arakawa_corner(b::CartesianBase{T}, i::Int, j::Int) where {T}
    (b.xlo + T(i - 1) * arakawa_dx(b),
     b.ylo + T(j - 1) * arakawa_dy(b))
end

_arakawa_base_esm(b::CartesianBase) = Dict{String,Any}(
    "family" => "cartesian",
    "nx" => b.nx,
    "ny" => b.ny,
    "extent" => [[b.xlo, b.ylo], [b.xhi, b.yhi]],
)

# ---------------------------------------------------------------------------
# ArakawaGrid: a base grid + stagger label
# ---------------------------------------------------------------------------

struct ArakawaGrid{T, B <: ArakawaBaseGrid} <: AbstractGrid
    base::B
    stagger::ArakawaStagger
    ghosts::Int
end

function ArakawaGrid(base::B, stagger::ArakawaStagger;
                     ghosts::Int = 0, dtype::Type{T} = Float64) where {T, B <: ArakawaBaseGrid}
    ghosts >= 0 || throw(DomainError(ghosts, "ghosts must be non-negative"))
    ArakawaGrid{T, B}(base, stagger, ghosts)
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
    arakawa_shape(g, loc)
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
    arakawa_cell_center(g.base, i, j)
end

"""
    u_face(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of the location where `u` lives at index `(i, j)` for this
stagger. Valid index range depends on the stagger (see `variable_shape`).
"""
function u_face(g::ArakawaGrid, i::Int, j::Int)
    _, u_loc, _ = arakawa_variable_locations(g.stagger)
    _check_bounds(g, u_loc, i, j)
    _coord_at(g, u_loc, i, j)
end

"""
    v_face(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of the location where `v` lives at index `(i, j)` for this
stagger.
"""
function v_face(g::ArakawaGrid, i::Int, j::Int)
    _, _, v_loc = arakawa_variable_locations(g.stagger)
    _check_bounds(g, v_loc, i, j)
    _coord_at(g, v_loc, i, j)
end

"""
    corners(g::ArakawaGrid, i::Int, j::Int) -> (x, y)

Coordinate of cell-corner `(i, j)` (`1 ≤ i ≤ nx+1`, `1 ≤ j ≤ ny+1`).
"""
function corners(g::ArakawaGrid, i::Int, j::Int)
    _check_bounds(g, Corner, i, j)
    arakawa_corner(g.base, i, j)
end

function _coord_at(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    loc === CellCenter ? arakawa_cell_center(g.base, i, j) :
    loc === UEdge      ? arakawa_x_edge(g.base, i, j) :
    loc === VEdge      ? arakawa_y_edge(g.base, i, j) :
    arakawa_corner(g.base, i, j)
end

function _check_bounds(g::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    ni, nj = arakawa_shape(g, loc)
    (1 <= i <= ni && 1 <= j <= nj) ||
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
    west  = i > 1  ? (loc, i - 1, j) : nothing
    east  = i < ni ? (loc, i + 1, j) : nothing
    south = j > 1  ? (loc, i, j - 1) : nothing
    north = j < nj ? (loc, i, j + 1) : nothing
    (west, east, south, north)
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
    if name === :dx
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
    Dict{String,Any}(
        "family"   => "arakawa",
        "dtype"    => _dtype_string(T),
        "topology" => "block_structured",
        "ghosts"   => g.ghosts,
        "n_cells"  => arakawa_nx(g) * arakawa_ny(g),
        "stagger"  => _stagger_string(g.stagger),
        "rotated"  => g.stagger === ArakawaE,
        "base"     => _arakawa_base_esm(g.base),
    )
end

function _stagger_string(s::ArakawaStagger)
    s === ArakawaA ? "A" :
    s === ArakawaB ? "B" :
    s === ArakawaC ? "C" :
    s === ArakawaD ? "D" : "E"
end

_dtype_string(::Type{Float64}) = "float64"
_dtype_string(::Type{Float32}) = "float32"
_dtype_string(::Type{T}) where {T} = string(T)
