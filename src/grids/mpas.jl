"""
MPAS unstructured Voronoi grid family (loader-backed).

Loader-backed per GRIDS_API.md §10 and mayor's 2026-04-20 scope correction on
dsc-3nw: `.esm` carries a small declarative config (family, dimensions,
loader ref). Geometry is derived on demand by the accessors
(`cell_centers`, `neighbors`, `cell_area`, `edge_length`, `metric_eval`).
Cross-binding conformance is compared at pinned query points, not serialized
bytes.

NetCDF I/O is intentionally not bundled — the caller must supply a
`reader_fn(path) -> MpasMeshData` when loading from a path. In-memory
construction via `mpas_mesh_data(...)` is the primary path for tests and
host-built meshes. Mirrors the Python `reader_fn` and TypeScript `readerFn`
contracts.

Index convention: Julia is the reference binding and uses 1-based cell and
edge indices with `0` as the "no neighbor" / "no edge" / "external boundary"
sentinel (per NetCDF convention). The Python / TypeScript / Rust bindings
use 0-based indices with `-1` sentinels; conformance compares accessor
outputs at pinned query points, not serialized adjacency arrays.
"""

struct MpasLoader
    path::String
    reader::String
    check::String
end

MpasLoader(;
    path::AbstractString, reader::AbstractString = "auto",
    check::AbstractString = "strict"
) =
    MpasLoader(String(path), String(reader), String(check))

const _MPAS_VALID_READERS = ("auto", "nc4", "mpas_mesh")
const _MPAS_VALID_CHECKS = ("strict", "lenient")

# Accept MpasLoader / NamedTuple / Dict as loader input.
_coerce_mpas_loader(ldr::MpasLoader) = _validate_mpas_loader(ldr)
function _coerce_mpas_loader(ldr::NamedTuple)
    path = get(ldr, :path, nothing)
    path === nothing && throw(ArgumentError("mpas: loader.path is required"))
    return _validate_mpas_loader(
        MpasLoader(;
            path = path,
            reader = String(get(ldr, :reader, "auto")),
            check = String(get(ldr, :check, "strict")),
        )
    )
end
function _coerce_mpas_loader(ldr::AbstractDict)
    path = get(ldr, "path", get(ldr, :path, nothing))
    path === nothing && throw(ArgumentError("mpas: loader.path is required"))
    return _validate_mpas_loader(
        MpasLoader(;
            path = String(path),
            reader = String(get(ldr, "reader", get(ldr, :reader, "auto"))),
            check = String(get(ldr, "check", get(ldr, :check, "strict"))),
        )
    )
end

function _validate_mpas_loader(ldr::MpasLoader)
    isempty(ldr.path) &&
        throw(ArgumentError("mpas: loader.path must be a non-empty string"))
    ldr.reader in _MPAS_VALID_READERS || throw(
        DomainError(
            ldr.reader,
            "mpas: loader.reader must be one of $(_MPAS_VALID_READERS)"
        )
    )
    ldr.check in _MPAS_VALID_CHECKS || throw(
        DomainError(
            ldr.check,
            "mpas: loader.check must be one of $(_MPAS_VALID_CHECKS)"
        )
    )
    return ldr
end

"""
    MpasMeshData

Validated in-memory MPAS Voronoi mesh (GRIDS_API.md §10, §7).

Field layout follows NetCDF / MPAS conventions. Julia stores the 2D
adjacency arrays column-major so the fast axis is the slot index
(`cells_on_cell[j, c]` = the `j`-th cell neighbor of cell `c`). Entries are
1-based cell/edge indices; `0` denotes "no neighbor" / "no edge" /
"external boundary".

Fields:
- `n_cells`, `n_edges`, `n_vertices`, `max_edges`: mesh dimensions
- `lon_cell`, `lat_cell`: (n_cells,) cell-center geographic coords (rad)
- `x_cell`, `y_cell`, `z_cell`: (n_cells,) cell-center cartesian (R-scaled)
- `area_cell`: (n_cells,) spherical-polygon area (m² if R in m)
- `n_edges_on_cell`: (n_cells,) valence of each cell (5 for pentagons, 6
  for hexagons in a quasi-uniform icosahedral-dual mesh)
- `cells_on_cell`: (max_edges, n_cells) 1-based neighbor cell indices
- `edges_on_cell`: (max_edges, n_cells) 1-based incident edge indices
- `lon_edge`, `lat_edge`: (n_edges,) edge-center geographic coords
- `cells_on_edge`: (2, n_edges) 1-based pair of cells sharing each edge
- `dc_edge`: (n_edges,) cell-center-to-cell-center great-circle distance
- `dv_edge`: (n_edges,) Voronoi-vertex-to-vertex great-circle arc length
"""
struct MpasMeshData
    n_cells::Int
    n_edges::Int
    n_vertices::Int
    max_edges::Int
    lon_cell::Vector{Float64}
    lat_cell::Vector{Float64}
    x_cell::Vector{Float64}
    y_cell::Vector{Float64}
    z_cell::Vector{Float64}
    area_cell::Vector{Float64}
    n_edges_on_cell::Vector{Int}
    cells_on_cell::Matrix{Int}
    edges_on_cell::Matrix{Int}
    lon_edge::Vector{Float64}
    lat_edge::Vector{Float64}
    cells_on_edge::Matrix{Int}
    dc_edge::Vector{Float64}
    dv_edge::Vector{Float64}
end

function _as_f64_1d(x, label::String, n::Int)
    v = collect(Float64, x)
    length(v) == n || throw(
        ArgumentError("mpas: $label must be length $n (got $(length(v)))")
    )
    return v
end

function _as_int_1d(x, label::String, n::Int)
    v = Vector{Int}(undef, n)
    i = 0
    for xi in x
        i += 1
        i > n && throw(ArgumentError("mpas: $label must be length $n (got > $n)"))
        v[i] = Int(xi)
    end
    i == n || throw(ArgumentError("mpas: $label must be length $n (got $i)"))
    return v
end

function _as_int_2d(x, label::String, rows::Int, cols::Int)
    if x isa AbstractMatrix
        size(x) == (rows, cols) || throw(
            ArgumentError(
                "mpas: $label must have shape ($rows, $cols); got $(size(x))"
            )
        )
        return Matrix{Int}(x)
    end
    # Accept flat iterable in column-major order.
    v = collect(x)
    length(v) == rows * cols || throw(
        ArgumentError(
            "mpas: $label must have $rows*$cols = $(rows * cols) entries " *
                "(got $(length(v)))"
        )
    )
    M = Matrix{Int}(undef, rows, cols)
    @inbounds for c in 1:cols, r in 1:rows
        M[r, c] = Int(v[(c - 1) * rows + r])
    end
    return M
end

"""
    mpas_mesh_data(; lon_cell, lat_cell, area_cell, n_edges_on_cell,
                     cells_on_cell, edges_on_cell,
                     lon_edge, lat_edge, cells_on_edge, dc_edge, dv_edge,
                     max_edges, x_cell=nothing, y_cell=nothing, z_cell=nothing,
                     n_vertices=0, R=6.371e6) -> MpasMeshData

Validated in-memory MPAS mesh constructor. Derives cartesian `x_cell`,
`y_cell`, `z_cell` from `(lon_cell, lat_cell, R)` when not supplied.

All adjacency arrays carry 1-based cell / edge indices; `0` is the "no
neighbor" / "no edge" / "external boundary" sentinel.
"""
function mpas_mesh_data(;
        lon_cell, lat_cell, area_cell, n_edges_on_cell,
        cells_on_cell, edges_on_cell,
        lon_edge, lat_edge, cells_on_edge, dc_edge, dv_edge,
        max_edges::Integer,
        x_cell = nothing, y_cell = nothing, z_cell = nothing,
        n_vertices::Integer = 0,
        R::Real = 6.371e6
    )
    max_edges = Int(max_edges)
    max_edges > 0 ||
        throw(DomainError(max_edges, "mpas: max_edges must be positive"))

    lon_cell_v = collect(Float64, lon_cell)
    n_cells = length(lon_cell_v)
    lon_edge_v = collect(Float64, lon_edge)
    n_edges = length(lon_edge_v)

    lat_cell_v = _as_f64_1d(lat_cell, "lat_cell", n_cells)
    area_cell_v = _as_f64_1d(area_cell, "area_cell", n_cells)
    neoc = _as_int_1d(n_edges_on_cell, "n_edges_on_cell", n_cells)
    coc = _as_int_2d(cells_on_cell, "cells_on_cell", max_edges, n_cells)
    eoc = _as_int_2d(edges_on_cell, "edges_on_cell", max_edges, n_cells)
    lat_edge_v = _as_f64_1d(lat_edge, "lat_edge", n_edges)
    coe = _as_int_2d(cells_on_edge, "cells_on_edge", 2, n_edges)
    dc_edge_v = _as_f64_1d(dc_edge, "dc_edge", n_edges)
    dv_edge_v = _as_f64_1d(dv_edge, "dv_edge", n_edges)

    R_f = Float64(R)
    x_cell_v = x_cell === nothing ?
        [R_f * cos(lat_cell_v[i]) * cos(lon_cell_v[i]) for i in 1:n_cells] :
        _as_f64_1d(x_cell, "x_cell", n_cells)
    y_cell_v = y_cell === nothing ?
        [R_f * cos(lat_cell_v[i]) * sin(lon_cell_v[i]) for i in 1:n_cells] :
        _as_f64_1d(y_cell, "y_cell", n_cells)
    z_cell_v = z_cell === nothing ?
        [R_f * sin(lat_cell_v[i]) for i in 1:n_cells] :
        _as_f64_1d(z_cell, "z_cell", n_cells)

    return MpasMeshData(
        n_cells, n_edges, Int(n_vertices), max_edges,
        lon_cell_v, lat_cell_v, x_cell_v, y_cell_v, z_cell_v, area_cell_v,
        neoc, coc, eoc,
        lon_edge_v, lat_edge_v, coe, dc_edge_v, dv_edge_v,
    )
end

"""
    check_mesh(mesh::MpasMeshData, strict::Bool)

Bounds-check adjacency arrays. In strict mode also enforces neighbor
reciprocity (if c lists n as a neighbor, n lists c).
"""
function check_mesh(mesh::MpasMeshData, strict::Bool)
    Nc = mesh.n_cells
    Ne = mesh.n_edges
    me = mesh.max_edges
    for c in 1:Nc
        k = mesh.n_edges_on_cell[c]
        (0 <= k <= me) || throw(
            AssertionError("mpas: n_edges_on_cell[$c]=$k out of [0, $me]")
        )
        for j in 1:k
            nb = mesh.cells_on_cell[j, c]
            (0 <= nb <= Nc) || throw(
                AssertionError(
                    "mpas: cells_on_cell[$j,$c]=$nb out of [0, $Nc]"
                )
            )
            e = mesh.edges_on_cell[j, c]
            (0 <= e <= Ne) || throw(
                AssertionError(
                    "mpas: edges_on_cell[$j,$c]=$e out of [0, $Ne]"
                )
            )
        end
    end
    for e in 1:Ne, s in 1:2
        c = mesh.cells_on_edge[s, e]
        (0 <= c <= Nc) || throw(
            AssertionError("mpas: cells_on_edge[$s,$e]=$c out of [0, $Nc]")
        )
    end
    strict || return nothing
    for c in 1:Nc
        k = mesh.n_edges_on_cell[c]
        for j in 1:k
            nb = mesh.cells_on_cell[j, c]
            nb == 0 && continue
            kb = mesh.n_edges_on_cell[nb]
            found = false
            for jj in 1:kb
                if mesh.cells_on_cell[jj, nb] == c
                    found = true
                    break
                end
            end
            found || throw(
                AssertionError(
                    "mpas: neighbor symmetry broken: cell $c -> $nb but " *
                        "not reverse"
                )
            )
        end
    end
    return nothing
end

"""
    MpasGrid{T} <: AbstractGrid

MPAS unstructured Voronoi grid. Wraps a validated `MpasMeshData` plus the
family options `(R, dtype, ghosts, loader)` and a provenance block (§6.4).
"""
struct MpasGrid{T} <: AbstractUnstructuredGrid
    mesh::MpasMeshData
    R::T
    dtype::String
    ghosts::Int
    loader::Union{MpasLoader, Nothing}
    provenance::Dict{String, Any}
end

const _MPAS_FAMILY_VERSION = "1.0.0"
const _MPAS_SOURCE_SHA = "dsc-7j0"
const _MPAS_READER_VERSION = "0.1.0"

"""
    build_mpas_grid(; mesh=nothing, loader=nothing, reader_fn=nothing,
                      R=6.371e6, dtype=Float64, ghosts=0) -> MpasGrid

Primary constructor for the MPAS family. One of `mesh` (an
`MpasMeshData`) or `loader` (an `MpasLoader` / NamedTuple / Dict with key
`path`) is required.

When loading from a path, `reader_fn(path::AbstractString) -> MpasMeshData`
must also be supplied — NetCDF I/O is not bundled with this package per
GRIDS_API.md §10.

`ghosts` must be 0 for loader-backed grids; the mesh file defines the full
interior.
"""
function build_mpas_grid(;
        mesh::Union{MpasMeshData, Nothing} = nothing,
        loader = nothing,
        reader_fn = nothing,
        R::Real = 6.371e6,
        dtype = Float64,
        ghosts::Int = 0
    )
    if !(dtype === Float64 || dtype === Float32)
        throw(ArgumentError("mpas: dtype must be Float64 or Float32, got $dtype"))
    end
    ghosts == 0 || throw(
        DomainError(
            ghosts,
            "mpas: ghosts must be 0 for loader-backed grids, got $ghosts"
        )
    )
    (R > 0 && isfinite(R)) ||
        throw(DomainError(R, "mpas: R must be a positive finite number"))

    local loader_record::Union{MpasLoader, Nothing}
    local resolved_mesh::MpasMeshData
    loader_record = nothing

    if mesh !== nothing
        resolved_mesh = mesh
        if loader !== nothing
            loader_record = _coerce_mpas_loader(loader)
        end
    elseif loader !== nothing
        loader_record = _coerce_mpas_loader(loader)
        reader_fn === nothing && throw(
            ArgumentError(
                "mpas: path-based loading requires reader_fn(path) -> " *
                    "MpasMeshData. NetCDF I/O is not bundled per " *
                    "GRIDS_API.md §10; pass reader_fn from the consumer."
            )
        )
        produced = reader_fn(loader_record.path)
        produced isa MpasMeshData || throw(
            ArgumentError("mpas: reader_fn must return an MpasMeshData")
        )
        resolved_mesh = produced
    else
        throw(
            ArgumentError(
                "mpas: provide either `mesh` (MpasMeshData) or `loader` " *
                    "with a `reader_fn`"
            )
        )
    end

    strict = loader_record === nothing ? true : loader_record.check == "strict"
    check_mesh(resolved_mesh, strict)

    T = dtype
    R_T = T(R)
    dtype_str = dtype === Float64 ? "float64" : "float32"

    provenance = Dict{String, Any}(
        "binding" => "julia",
        "family" => "mpas",
        "version" => _MPAS_FAMILY_VERSION,
        "source_sha" => _MPAS_SOURCE_SHA,
        "reader_version" => _MPAS_READER_VERSION,
        "dtype" => dtype_str,
    )
    if loader_record !== nothing
        provenance["loader"] = Dict{String, Any}(
            "path" => loader_record.path,
            "reader" => loader_record.reader,
            "check" => loader_record.check,
        )
    end

    return MpasGrid{T}(
        resolved_mesh, R_T, dtype_str, ghosts, loader_record, provenance
    )
end

# --- Accessors ---------------------------------------------------------

"""
    cell_centers(grid::MpasGrid) -> NamedTuple{(:lon,:lat)}
    cell_centers(grid::MpasGrid, c::Integer) -> (lon, lat)
"""
cell_centers(g::MpasGrid) = (lon = g.mesh.lon_cell, lat = g.mesh.lat_cell)
function cell_centers(g::MpasGrid, c::Integer)
    _check_mpas_cell(g, c)
    return (g.mesh.lon_cell[c], g.mesh.lat_cell[c])
end

"""
    cell_center_cart(grid::MpasGrid, c::Integer) -> (x, y, z)

Cartesian cell-center coords (R-scaled).
"""
function cell_center_cart(g::MpasGrid, c::Integer)
    _check_mpas_cell(g, c)
    return (g.mesh.x_cell[c], g.mesh.y_cell[c], g.mesh.z_cell[c])
end

"""
    neighbors(grid::MpasGrid, c::Integer) -> Vector{Int}

Returns 1-based cell indices sharing an edge with cell `c`. Boundary
sentinels (`0`) are filtered out; closed global meshes return a k-tuple
where k is the cell's valence (5 for pentagons, 6 for hexagons).
"""
function neighbors(g::MpasGrid, c::Integer)
    _check_mpas_cell(g, c)
    k = g.mesh.n_edges_on_cell[c]
    out = Int[]
    sizehint!(out, k)
    @inbounds for j in 1:k
        nb = g.mesh.cells_on_cell[j, c]
        nb != 0 && push!(out, nb)
    end
    return out
end

"""
    cell_area(grid::MpasGrid, c::Integer) -> Float64
"""
function cell_area(g::MpasGrid, c::Integer)
    _check_mpas_cell(g, c)
    return g.mesh.area_cell[c]
end

"""
    edge_length(grid::MpasGrid, e::Integer) -> Float64

Voronoi-vertex-to-Voronoi-vertex arc length (`dv_edge`).
"""
function edge_length(g::MpasGrid, e::Integer)
    _check_mpas_edge(g, e)
    return g.mesh.dv_edge[e]
end

"""
    cell_distance(grid::MpasGrid, e::Integer) -> Float64

Cell-center-to-cell-center great-circle distance (`dc_edge`).
"""
function cell_distance(g::MpasGrid, e::Integer)
    _check_mpas_edge(g, e)
    return g.mesh.dc_edge[e]
end

"""
    metric_eval(grid::MpasGrid, name::Symbol, i::Integer)

Per-cell metrics: `:lon`, `:lat`, `:area`, `:x`, `:y`, `:z`,
`:n_edges_on_cell`. Per-edge metrics: `:lon_edge`, `:lat_edge`,
`:dc_edge`, `:dv_edge`.
"""
function metric_eval(g::MpasGrid, name::Symbol, i::Integer)
    if name === :lon
        _check_mpas_cell(g, i); return g.mesh.lon_cell[i]
    elseif name === :lat
        _check_mpas_cell(g, i); return g.mesh.lat_cell[i]
    elseif name === :area
        _check_mpas_cell(g, i); return g.mesh.area_cell[i]
    elseif name === :x
        _check_mpas_cell(g, i); return g.mesh.x_cell[i]
    elseif name === :y
        _check_mpas_cell(g, i); return g.mesh.y_cell[i]
    elseif name === :z
        _check_mpas_cell(g, i); return g.mesh.z_cell[i]
    elseif name === :n_edges_on_cell
        _check_mpas_cell(g, i); return Float64(g.mesh.n_edges_on_cell[i])
    elseif name === :lon_edge
        _check_mpas_edge(g, i); return g.mesh.lon_edge[i]
    elseif name === :lat_edge
        _check_mpas_edge(g, i); return g.mesh.lat_edge[i]
    elseif name === :dc_edge
        _check_mpas_edge(g, i); return g.mesh.dc_edge[i]
    elseif name === :dv_edge
        _check_mpas_edge(g, i); return g.mesh.dv_edge[i]
    else
        throw(ArgumentError("mpas: metric_eval: unknown metric $name"))
    end
end

metric_eval(g::MpasGrid, name::AbstractString, i::Integer) =
    metric_eval(g, Symbol(name), i)

n_cells(g::MpasGrid) = g.mesh.n_cells
n_edges(g::MpasGrid) = g.mesh.n_edges
n_vertices(g::MpasGrid) = g.mesh.n_vertices
max_edges(g::MpasGrid) = g.mesh.max_edges
total_area(g::MpasGrid) = sum(g.mesh.area_cell)
family(::MpasGrid) = "mpas"

function _check_mpas_cell(g::MpasGrid, c::Integer)
    (1 <= c <= g.mesh.n_cells) || throw(
        BoundsError(g.mesh.lon_cell, c)
    )
end

function _check_mpas_edge(g::MpasGrid, e::Integer)
    (1 <= e <= g.mesh.n_edges) || throw(
        BoundsError(g.mesh.lon_edge, e)
    )
end

# --- .esm lowering (small declarative config per mayor's correction) ---

"""
    to_esm(grid::MpasGrid) -> Dict{String,Any}

Returns a §6-schema-valid declarative config: family + dimensions +
loader ref + provenance. No inline geometry arrays.
"""
function to_esm(g::MpasGrid)
    loader_blk = g.loader === nothing ? nothing :
        Dict{String, Any}(
            "path" => g.loader.path,
            "reader" => g.loader.reader,
            "check" => g.loader.check,
        )
    return Dict{String, Any}(
        "family" => "mpas",
        "version" => _MPAS_FAMILY_VERSION,
        "schema_version" => _MPAS_FAMILY_VERSION,
        "topology" => "unstructured",
        "dtype" => g.dtype,
        "ghosts" => g.ghosts,
        "n_cells" => n_cells(g),
        "n_edges" => n_edges(g),
        "n_vertices" => n_vertices(g),
        "max_edges" => max_edges(g),
        "options" => Dict{String, Any}(
            "R" => Float64(g.R),
            "loader" => loader_blk,
        ),
        "provenance" => g.provenance,
    )
end

# ---------------------------------------------------------------------------
# ESS Grid trait — Tier C + Tier U (unstructured Voronoi sphere)
#
# MPAS uses a flat 1-D cell axis (`:cell`); spatial coordinates are exposed
# via `cell_centers(g, :lon)` / `(g, :lat)` (still flat-vector, length
# `n_cells`). Per-cell neighbour slots are reached through
# `neighbor_indices(g, :cell, k)` for `k ∈ 1:max_edges` — the offset is the
# slot index, not an axial step. Boundary sentinel: `0`. Tier-U methods
# (`cell_neighbor_table`, `cell_valence`, `edge_length`, `cell_distance`)
# expose the ragged adjacency directly.
# ---------------------------------------------------------------------------

n_dims(::MpasGrid) = 1
axis_names(::MpasGrid) = (:cell,)

function _mpas_axis_check(::MpasGrid, axis::Symbol)
    axis in (:cell, :lon, :lat) ||
        throw(ArgumentError("mpas: unknown axis :$axis (expected :cell, :lon, or :lat)"))
    return nothing
end

function cell_centers(g::MpasGrid{T}, axis::Symbol) where {T}
    _mpas_axis_check(g, axis)
    axis === :cell &&
        throw(ArgumentError("mpas: cell_centers needs a coordinate axis (:lon or :lat); :cell is the layout axis"))
    return _grid_memo!(g, (:cell_centers, axis)) do
        src = axis === :lon ? g.mesh.lon_cell : g.mesh.lat_cell
        return T[T(x) for x in src]
    end
end

function cell_widths(g::MpasGrid{T}, axis::Symbol) where {T}
    _mpas_axis_check(g, axis)
    # MPAS Voronoi cells are non-axial. Expose an isotropic length proxy:
    # `sqrt(area_cell)`. Consumers needing edge lengths or center-to-center
    # distances should use `edge_length(g, e)` / `cell_distance(g, e)`.
    return _grid_memo!(g, (:cell_widths, axis)) do
        out = Vector{T}(undef, g.mesh.n_cells)
        @inbounds for c in 1:g.mesh.n_cells
            out[c] = T(sqrt(g.mesh.area_cell[c]))
        end
        return out
    end
end

function cell_volume(g::MpasGrid{T}) where {T}
    return _grid_memo!(g, :cell_volume) do
        T[T(x) for x in g.mesh.area_cell]
    end
end

function neighbor_indices(g::MpasGrid, axis::Symbol, offset::Int)
    axis === :cell ||
        throw(ArgumentError("mpas: neighbor_indices: only axis=:cell is supported (offset = neighbour slot 1..max_edges); got :$axis"))
    me = g.mesh.max_edges
    (1 <= offset <= me) ||
        throw(ArgumentError("mpas: neighbor_indices offset (slot index) must lie in 1..max_edges=$me; got $offset"))
    return _grid_memo!(g, (:neighbor_indices, axis, offset)) do
        Nc = g.mesh.n_cells
        out = Vector{Int}(undef, Nc)
        @inbounds for c in 1:Nc
            valence = g.mesh.n_edges_on_cell[c]
            out[c] = offset <= valence ? g.mesh.cells_on_cell[offset, c] : 0
        end
        return out
    end
end

function boundary_mask(g::MpasGrid, axis::Symbol, side::Symbol)
    axis === :cell ||
        throw(ArgumentError("mpas: boundary_mask: only axis=:cell is supported; got :$axis"))
    side in (:lower, :upper) ||
        throw(ArgumentError("mpas: side must be :lower or :upper; got :$side"))
    # A "boundary" cell on an MPAS mesh is any cell with at least one zero
    # entry in `cells_on_cell` — i.e. an open-domain edge. Closed global
    # meshes have none. The `side` argument is preserved for trait shape but
    # both sides return the same flag (no axial direction to distinguish).
    return _grid_memo!(g, (:boundary_mask, axis, side)) do
        Nc = g.mesh.n_cells
        out = falses(Nc)
        @inbounds for c in 1:Nc
            v = g.mesh.n_edges_on_cell[c]
            for j in 1:v
                if g.mesh.cells_on_cell[j, c] == 0
                    out[c] = true
                    break
                end
            end
        end
        return out
    end
end

# Tier U — ragged adjacency.

"""
    cell_neighbor_table(g::MpasGrid) -> Matrix{Int}

`(max_edges, n_cells)` neighbour-cell table. Slot index is the local edge
index on each cell; `0` denotes "no neighbour" (boundary cell on open
meshes).
"""
cell_neighbor_table(g::MpasGrid) = g.mesh.cells_on_cell

"""
    cell_valence(g::MpasGrid) -> Vector{Int}

Per-cell number of neighbouring cells (`n_edges_on_cell`).
"""
cell_valence(g::MpasGrid) = g.mesh.n_edges_on_cell

# Tier-U `edge_length` and `cell_distance` already exist as scalar accessors
# `(g, e::Integer)`. Provide bulk-array siblings that return the full edge
# vectors per the trait contract.

"""
    edge_length(g::MpasGrid) -> Vector{Float64}

Voronoi-vertex-to-vertex arc lengths (`dv_edge`), shape `(n_edges,)`.
"""
edge_length(g::MpasGrid) = g.mesh.dv_edge

"""
    cell_distance(g::MpasGrid) -> Vector{Float64}

Cell-center-to-cell-center great-circle distances (`dc_edge`), shape
`(n_edges,)`.
"""
cell_distance(g::MpasGrid) = g.mesh.dc_edge
