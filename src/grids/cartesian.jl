"""
Cartesian grid family (1D / 2D / 3D, uniform + non-uniform).

Implements the public API contract defined in `docs/GRIDS_API.md`. The
generator function lives in the `EarthSciDiscretizations.grids` submodule
(see bottom of this file) and is invoked as:

    using EarthSciDiscretizations
    g = EarthSciDiscretizations.grids.cartesian(; nx=10, extent=[(0.0, 1.0)])
    g = EarthSciDiscretizations.grids.cartesian(; nx=8, ny=4,
                                                  extent=[(0.0, 1.0), (0.0, 0.5)])
    g = EarthSciDiscretizations.grids.cartesian(; edges=[xe, ye])     # non-uniform

Per the bead correction (mayor 2026-04-20): the on-disk `.esm` form is a
small declarative config (family + dimensions + extents), NOT a serialized
geometry blob. This runtime provides accessors (`cell_centers`,
`cell_widths`, `cell_volume`, `neighbors`, `metric_eval`) that DERIVE
geometry from that declaration via pure math.
"""

# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

struct CartesianGrid{T <: AbstractFloat, N} <: AbstractCurvilinearGrid
    n::NTuple{N, Int}
    extent::NTuple{N, Tuple{T, T}}
    edges::NTuple{N, Vector{T}}
    centers::NTuple{N, Vector{T}}
    widths::NTuple{N, Vector{T}}
    uniform::NTuple{N, Bool}
    ghosts::Int
end

Base.ndims(::CartesianGrid{T, N}) where {T, N} = N
Base.eltype(::CartesianGrid{T, N}) where {T, N} = T

# ---------------------------------------------------------------------------
# Constructor helpers
# ---------------------------------------------------------------------------

const _AXIS_NAMES = (:x, :y, :z)

function _normalize_extent(extent, ::Type{T}, ndim::Int) where {T}
    if extent isa AbstractMatrix
        size(extent, 1) == 2 ||
            throw(ArgumentError("cartesian: extent matrix must have shape (2, ndim); got $(size(extent))"))
        size(extent, 2) == ndim ||
            throw(ArgumentError("cartesian: extent matrix has $(size(extent, 2)) cols, expected $(ndim)"))
        return ntuple(d -> (T(extent[1, d]), T(extent[2, d])), ndim)
    elseif extent isa AbstractVector
        length(extent) == ndim ||
            throw(ArgumentError("cartesian: extent vector has length $(length(extent)), expected $(ndim)"))
        out = ntuple(ndim) do d
            e = extent[d]
            (e isa Tuple || e isa AbstractVector) && length(e) == 2 ||
                throw(ArgumentError("cartesian: extent[$d] must be a 2-element (lo, hi); got $e"))
            (T(e[1]), T(e[2]))
        end
        return out
    else
        throw(ArgumentError("cartesian: extent must be a Matrix or Vector of (lo, hi); got $(typeof(extent))"))
    end
end

function _uniform_axis(n::Int, lo::T, hi::T) where {T}
    n ≥ 1 || throw(ArgumentError("cartesian: cell count per axis must be ≥ 1; got $n"))
    hi > lo || throw(DomainError((lo, hi), "cartesian: extent must satisfy hi > lo"))
    dx = (hi - lo) / n
    edges = T[lo + i * dx for i in 0:n]
    centers = T[(edges[i] + edges[i + 1]) / 2 for i in 1:n]
    widths = T[edges[i + 1] - edges[i] for i in 1:n]
    return edges, centers, widths
end

function _nonuniform_axis(edges_in::AbstractVector, ::Type{T}) where {T}
    length(edges_in) ≥ 2 ||
        throw(ArgumentError("cartesian: non-uniform edges must have ≥ 2 entries; got $(length(edges_in))"))
    edges = T[T(e) for e in edges_in]
    for i in 1:(length(edges) - 1)
        edges[i + 1] > edges[i] ||
            throw(DomainError(edges, "cartesian: edges must be strictly increasing"))
    end
    n = length(edges) - 1
    centers = T[(edges[i] + edges[i + 1]) / 2 for i in 1:n]
    widths = T[edges[i + 1] - edges[i] for i in 1:n]
    return n, edges, centers, widths
end

function _is_uniform(widths::Vector{T}) where {T}
    length(widths) ≤ 1 && return true
    w0 = widths[1]
    tol = max(eps(T), eps(T) * abs(w0)) * 8
    return all(w -> abs(w - w0) ≤ tol, widths)
end

# ---------------------------------------------------------------------------
# Generator (entry point)
# ---------------------------------------------------------------------------

function _cartesian(;
        nx::Union{Int, Nothing} = nothing,
        ny::Union{Int, Nothing} = nothing,
        nz::Union{Int, Nothing} = nothing,
        extent = nothing,
        edges = nothing,
        dtype::Type{T} = Float64,
        ghosts::Int = 0
    ) where {T <: AbstractFloat}

    ghosts ≥ 0 || throw(ArgumentError("cartesian: ghosts must be ≥ 0; got $ghosts"))

    if edges !== nothing
        # Non-uniform path.
        (extent === nothing && nx === nothing && ny === nothing && nz === nothing) ||
            throw(ArgumentError("cartesian: pass either `edges` (non-uniform) OR `nx/ny/nz`+`extent` (uniform), not both"))
        edges isa AbstractVector && !isempty(edges) ||
            throw(ArgumentError("cartesian: `edges` must be a non-empty Vector of edge arrays"))
        ndim = length(edges)
        ndim ≤ 3 || throw(ArgumentError("cartesian: ndim must be ≤ 3; got $ndim"))
        ax_n = Vector{Int}(undef, ndim)
        ax_e = Vector{Vector{T}}(undef, ndim)
        ax_c = Vector{Vector{T}}(undef, ndim)
        ax_w = Vector{Vector{T}}(undef, ndim)
        ax_u = Vector{Bool}(undef, ndim)
        ax_x = Vector{Tuple{T, T}}(undef, ndim)
        for d in 1:ndim
            n_d, e_d, c_d, w_d = _nonuniform_axis(edges[d], T)
            ax_n[d] = n_d
            ax_e[d] = e_d
            ax_c[d] = c_d
            ax_w[d] = w_d
            ax_u[d] = _is_uniform(w_d)
            ax_x[d] = (e_d[1], e_d[end])
        end
        return CartesianGrid{T, ndim}(
            ntuple(d -> ax_n[d], ndim),
            ntuple(d -> ax_x[d], ndim),
            ntuple(d -> ax_e[d], ndim),
            ntuple(d -> ax_c[d], ndim),
            ntuple(d -> ax_w[d], ndim),
            ntuple(d -> ax_u[d], ndim),
            ghosts,
        )
    end

    # Uniform path: nx (and optionally ny, nz) + extent.
    nx === nothing && throw(ArgumentError("cartesian: required option `nx` is missing"))
    ns = (nx, ny, nz)
    ndim = ny === nothing ? 1 : (nz === nothing ? 2 : 3)
    if ndim ≥ 2
        ny === nothing && throw(ArgumentError("cartesian: required option `ny` is missing for 2D/3D grid"))
    end
    if ndim == 3
        nz === nothing && throw(ArgumentError("cartesian: required option `nz` is missing for 3D grid"))
    end
    extent === nothing && throw(ArgumentError("cartesian: required option `extent` is missing"))
    ext_tuple = _normalize_extent(extent, T, ndim)

    ax_e = Vector{Vector{T}}(undef, ndim)
    ax_c = Vector{Vector{T}}(undef, ndim)
    ax_w = Vector{Vector{T}}(undef, ndim)
    for d in 1:ndim
        e, c, w = _uniform_axis(ns[d], ext_tuple[d][1], ext_tuple[d][2])
        ax_e[d] = e; ax_c[d] = c; ax_w[d] = w
    end

    return CartesianGrid{T, ndim}(
        ntuple(d -> ns[d], ndim),
        ext_tuple,
        ntuple(d -> ax_e[d], ndim),
        ntuple(d -> ax_c[d], ndim),
        ntuple(d -> ax_w[d], ndim),
        ntuple(_ -> true, ndim),
        ghosts,
    )
end

# ---------------------------------------------------------------------------
# Accessors (the runtime contract per the mayor's correction)
# ---------------------------------------------------------------------------

"""
    cell_centers(grid::CartesianGrid, i, [j, [k]])

Return the geometric center of the cell at the given index as an
`NTuple{N,T}`. Pure math from extents — no stored geometry blob.
"""
function cell_centers(grid::CartesianGrid{T, N}, idx::Vararg{Int, N}) where {T, N}
    return ntuple(d -> grid.centers[d][idx[d]], N)
end

"""
    cell_widths(grid::CartesianGrid, i, [j, [k]])

Return the per-axis cell widths at the given index as an `NTuple{N,T}`.
"""
function cell_widths(grid::CartesianGrid{T, N}, idx::Vararg{Int, N}) where {T, N}
    return ntuple(d -> grid.widths[d][idx[d]], N)
end

"""
    cell_volume(grid::CartesianGrid, i, [j, [k]])

Return the cell measure: length (1D), area (2D), or volume (3D).
"""
function cell_volume(grid::CartesianGrid{T, N}, idx::Vararg{Int, N}) where {T, N}
    return prod(cell_widths(grid, idx...))
end

"""
    neighbors(grid::CartesianGrid, i, [j, [k]])

Return the axis-aligned neighbors of a cell as a Dict mapping
`(axis, side)` → neighbor index tuple. `axis ∈ 1:N`, `side ∈ (-1, +1)`.
A boundary cell with no neighbor on a given side is omitted from the
returned mapping (rather than carrying a `nothing`), keeping the
return type concrete.
"""
function neighbors(grid::CartesianGrid{T, N}, idx::Vararg{Int, N}) where {T, N}
    out = Dict{Tuple{Int, Int}, NTuple{N, Int}}()
    for d in 1:N
        i_d = idx[d]
        if i_d > 1
            out[(d, -1)] = ntuple(k -> k == d ? i_d - 1 : idx[k], N)
        end
        if i_d < grid.n[d]
            out[(d, +1)] = ntuple(k -> k == d ? i_d + 1 : idx[k], N)
        end
    end
    return out
end

"""
    metric_eval(grid::CartesianGrid, name::Symbol, i, [j, [k]])

Evaluate a named metric at a cell. Defined names:

- `:volume` — cell measure (length / area / volume).
- `:dx`, `:dy`, `:dz` — per-axis widths (1D / 2D / 3D as applicable).
- `:face_area_x` (`:_y`, `:_z`) — face area normal to the named axis.
  In 1D the face area is `1`; in 2D it is the orthogonal cell width;
  in 3D it is the product of the two orthogonal widths.
- `:g` — metric tensor (identity for cartesian) as an `NTuple{N,NTuple{N,T}}`.
- `:jacobian` — `det(g)^{1/2}` = 1 for cartesian, returned as `T(1)`.
"""
function metric_eval(grid::CartesianGrid{T, N}, name::Symbol, idx::Vararg{Int, N}) where {T, N}
    name === :volume && return cell_volume(grid, idx...)
    name === :jacobian && return one(T)
    if name === :g
        return ntuple(i -> ntuple(j -> i == j ? one(T) : zero(T), N), N)
    end
    # axis widths
    if name === :dx
        _check_axis(N, 1); return grid.widths[1][idx[1]]
    end
    if name === :dy
        _check_axis(N, 2); return grid.widths[2][idx[2]]
    end
    if name === :dz
        _check_axis(N, 3); return grid.widths[3][idx[3]]
    end
    if name === :face_area_x
        return _face_area(grid, 1, idx)
    end
    if name === :face_area_y
        _check_axis(N, 2); return _face_area(grid, 2, idx)
    end
    if name === :face_area_z
        _check_axis(N, 3); return _face_area(grid, 3, idx)
    end
    throw(ArgumentError("cartesian: unknown metric name :$name"))
end

_check_axis(N::Int, d::Int) =
    d ≤ N || throw(ArgumentError("cartesian: axis $d not present in $(N)D grid"))

function _face_area(grid::CartesianGrid{T, N}, axis::Int, idx::NTuple{N, Int}) where {T, N}
    # Product of widths on all axes EXCEPT `axis`. In 1D this is empty → 1.
    a = one(T)
    for d in 1:N
        d == axis && continue
        a *= grid.widths[d][idx[d]]
    end
    return a
end

# ---------------------------------------------------------------------------
# Provenance + .esm lowering
# ---------------------------------------------------------------------------

const _BINDING_VERSION = "1.0.0"
const _API_VERSION = "1.0.0"

# `_dtype_string` is shared with other grid families (e.g., arakawa) at
# the parent module scope; do NOT redefine here to avoid method-overwriting
# warnings during precompile.

function _math_lib_fingerprint()
    if Sys.islinux()
        return "libm-glibc"
    elseif Sys.isapple()
        return "libm-apple"
    elseif Sys.iswindows()
        return "libm-msvc"
    else
        return "libm-unknown"
    end
end

function _provenance(::CartesianGrid)
    return Dict{String, Any}(
        "binding" => "julia",
        "binding_version" => _BINDING_VERSION,
        "api_version" => _API_VERSION,
        "platform" => string(Sys.MACHINE),
        "runtime" => "julia-" * string(VERSION),
        "math_lib" => _math_lib_fingerprint(),
        "source_sha" => "",
    )
end

"""
    to_esm(grid::CartesianGrid) -> Dict{String,Any}

Lower a cartesian grid to the §6-schema declarative form. Per the mayor's
correction (2026-04-20), this is a small config (family / dims / extents /
optionally explicit edges for non-uniform), NOT a serialized geometry
payload. Cross-binding conformance compares accessor outputs at pinned
query points, with this `.esm` payload as the round-trip ground truth.
"""
function to_esm(grid::CartesianGrid{T, N}) where {T, N}
    d = Dict{String, Any}(
        "family" => "cartesian",
        "version" => _API_VERSION,
        "dtype" => _dtype_string(T),
        "topology" => "rectilinear",
        "ndim" => N,
        "ghosts" => grid.ghosts,
        "n_cells" => prod(grid.n),
        "n" => collect(grid.n),
        "extent" => Vector{Any}([Any[grid.extent[a][1], grid.extent[a][2]] for a in 1:N]),
        "uniform" => collect(grid.uniform),
        "provenance" => _provenance(grid),
    )
    # Non-uniform axes need their explicit edges to round-trip.
    if any(.!grid.uniform)
        d["edges"] = Vector{Any}(
            [
                grid.uniform[a] ? Any[] : Any[Float64(e) for e in grid.edges[a]]
                    for a in 1:N
            ]
        )
    end
    return d
end

# ---------------------------------------------------------------------------
# Ecosystem integration stub (Meshes.jl primary — optional)
# ---------------------------------------------------------------------------

"""
    as_meshes(grid::CartesianGrid)

Return a `Meshes.CartesianGrid` view of the grid. Requires `Meshes.jl` to
be loaded; if not available, an `ArgumentError` instructs the caller. Per
GRIDS_API.md §5.1, the fallback (no Meshes.jl) is always usable through
the accessors above.
"""
function as_meshes(::CartesianGrid)
    throw(
        ArgumentError(
            "as_meshes requires Meshes.jl. Install and `using Meshes` to enable; " *
                "the cartesian generator + accessors work without it."
        )
    )
end

# Public namespace registration: the `cartesian` entry point in the
# `EarthSciDiscretizations.grids` submodule (defined at the bottom of
# `src/EarthSciDiscretizations.jl`) aliases `_cartesian` here. See
# GRIDS_API.md §2.3 for the public call form.

# ---------------------------------------------------------------------------
# ESS Grid trait — Tier C + Tier M (identity metric)
#
# Cartesian grids subtype `AbstractCurvilinearGrid` and supply an identity
# metric per RFC §1 Tier M, so they can drop into curvilinear assemblers
# (`precompute_laplacian_stencil`, `precompute_gradient_stencil`) unchanged.
# ---------------------------------------------------------------------------

n_dims(::CartesianGrid{T, N}) where {T, N} = N
n_cells(g::CartesianGrid) = prod(g.n)
axis_names(::CartesianGrid{T, N}) where {T, N} = ntuple(d -> _AXIS_NAMES[d], N)

function _cartesian_axis_idx(::CartesianGrid{T, N}, axis::Symbol) where {T, N}
    @inbounds for d in 1:N
        _AXIS_NAMES[d] === axis && return d
    end
    throw(ArgumentError("cartesian: unknown axis :$axis for $(N)D grid"))
end

function cell_centers(g::CartesianGrid{T, N}, axis::Symbol) where {T, N}
    d = _cartesian_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_centers, axis)) do
        nc = prod(g.n)
        out = Vector{T}(undef, nc)
        cs = g.centers[d]
        ci = CartesianIndices(g.n)
        @inbounds for k in 1:nc
            out[k] = cs[ci[k][d]]
        end
        return out
    end
end

function cell_widths(g::CartesianGrid{T, N}, axis::Symbol) where {T, N}
    d = _cartesian_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_widths, axis)) do
        nc = prod(g.n)
        out = Vector{T}(undef, nc)
        ws = g.widths[d]
        ci = CartesianIndices(g.n)
        @inbounds for k in 1:nc
            out[k] = ws[ci[k][d]]
        end
        return out
    end
end

function cell_volume(g::CartesianGrid{T, N}) where {T, N}
    return _grid_memo!(g, :cell_volume) do
        nc = prod(g.n)
        out = Vector{T}(undef, nc)
        ci = CartesianIndices(g.n)
        @inbounds for k in 1:nc
            v = one(T)
            ix = ci[k]
            for d in 1:N
                v *= g.widths[d][ix[d]]
            end
            out[k] = v
        end
        return out
    end
end

function neighbor_indices(g::CartesianGrid{T, N}, axis::Symbol, offset::Int) where {T, N}
    d = _cartesian_axis_idx(g, axis)
    return _grid_memo!(g, (:neighbor_indices, axis, offset)) do
        nc = prod(g.n)
        out = Vector{Int}(undef, nc)
        ci = CartesianIndices(g.n)
        li = LinearIndices(g.n)
        nd = g.n[d]
        @inbounds for k in 1:nc
            ix = ci[k]
            new_d = ix[d] + offset
            if 1 <= new_d <= nd
                new_ix = ntuple(p -> p == d ? new_d : ix[p], N)
                out[k] = li[CartesianIndex(new_ix)]
            else
                out[k] = 0
            end
        end
        return out
    end
end

function boundary_mask(g::CartesianGrid{T, N}, axis::Symbol, side::Symbol) where {T, N}
    d = _cartesian_axis_idx(g, axis)
    side in (:lower, :upper) ||
        throw(ArgumentError("cartesian: side must be :lower or :upper; got :$side"))
    return _grid_memo!(g, (:boundary_mask, axis, side)) do
        nc = prod(g.n)
        out = Vector{Bool}(undef, nc)
        target = side === :lower ? 1 : g.n[d]
        ci = CartesianIndices(g.n)
        @inbounds for k in 1:nc
            out[k] = ci[k][d] == target
        end
        return out
    end
end

# Tier-M identity metric — `g_ij = δ_ij`, `J = ∏ dx_d`, derivatives of `g_ij`
# vanish. The coordinate Jacobian is identity vs the cartesian target, and
# the second-derivative Jacobian is zero.

function metric_g(g::CartesianGrid{T, N}) where {T, N}
    return _grid_memo!(g, :metric_g) do
        nc = prod(g.n)
        out = zeros(T, nc, N, N)
        @inbounds for k in 1:nc, d in 1:N
            out[k, d, d] = one(T)
        end
        return out
    end
end

function metric_ginv(g::CartesianGrid{T, N}) where {T, N}
    return _grid_memo!(g, :metric_ginv) do
        nc = prod(g.n)
        out = zeros(T, nc, N, N)
        @inbounds for k in 1:nc, d in 1:N
            out[k, d, d] = one(T)
        end
        return out
    end
end

function metric_jacobian(g::CartesianGrid{T, N}) where {T, N}
    # J = product of cell widths (cell_volume on a cartesian grid).
    return cell_volume(g)
end

function metric_dgij_dxk(g::CartesianGrid{T, N}) where {T, N}
    return _grid_memo!(g, :metric_dgij_dxk) do
        zeros(T, prod(g.n), N, N, N)
    end
end

function coord_jacobian(g::CartesianGrid{T, N}, target::Symbol) where {T, N}
    target === :cartesian ||
        throw(ArgumentError("cartesian: coord_jacobian only supports target=:cartesian; got :$target"))
    return _grid_memo!(g, (:coord_jacobian, target)) do
        nc = prod(g.n)
        out = zeros(T, nc, N, N)
        @inbounds for k in 1:nc, d in 1:N
            out[k, d, d] = one(T)
        end
        return out
    end
end

function coord_jacobian_second(g::CartesianGrid{T, N}, target::Symbol) where {T, N}
    target === :cartesian ||
        throw(ArgumentError("cartesian: coord_jacobian_second only supports target=:cartesian; got :$target"))
    return _grid_memo!(g, (:coord_jacobian_second, target)) do
        zeros(T, prod(g.n), N, N, N)
    end
end
