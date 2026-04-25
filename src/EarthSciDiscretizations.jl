module EarthSciDiscretizations

using LinearAlgebra: norm, cross, dot
using SymbolicUtils: SymReal, BSImpl, idxs_for_arrayop, @arrayop, @syms
using SymbolicUtils
using Symbolics: unwrap, wrap, Differential, iscall, operation, arguments, Num
import Symbolics

# Grid infrastructure
include("grids/abstract_grid.jl")
include("grids/panel_connectivity.jl")
include("grids/metric_tensors.jl")
include("grids/super_grid.jl")
include("grids/cubed_sphere.jl")
include("grids/duo.jl")
include("grids/mpas.jl")
include("grids/cartesian.jl")
include("grids/vertical.jl")
include("grids/latlon.jl")

# Staggering and discrete space
include("staggering.jl")
include("discrete_space.jl")
include("ghost_cells.jl")

# Arakawa staggering runtime (depends on VarLocation from staggering.jl)
include("grids/arakawa.jl")

# Operator utilities
include("operators/arrayop_utils.jl")

# FV operators (all ArrayOp-based)
include("operators/divergence.jl")
include("operators/gradient.jl")
include("operators/laplacian.jl")
include("operators/reconstruction.jl")
include("operators/flux_1d.jl")
include("operators/transport_2d.jl")
include("operators/vertical_remap.jl")

# FV3-specific operators
include("operators/wind_ops.jl")
include("operators/vorticity.jl")
include("operators/kinetic_energy.jl")
include("operators/ppm_edge.jl")

# Boundary-condition handler
include("bc_handler.jl")

# Discretization rule catalog (parser delegator; ESS integration pending)
include("rules.jl")

# AST coefficient evaluator (thin passthrough to EarthSciSerialization)
include("rule_eval.jl")

# Exports: Grid types
export AbstractGrid, AbstractCubedSphereGrid, CubedSphereGrid
export AbstractCurvilinearGrid, AbstractStaggeredGrid,
    AbstractVerticalGrid, AbstractUnstructuredGrid
export GridTraitError
export total_area

# Exports: ESS Grid trait (Tier C / Tier M)
export n_dims, axis_names, neighbor_indices, boundary_mask
export metric_g, metric_ginv, metric_jacobian, metric_dgij_dxk
export coord_jacobian, coord_jacobian_second

# Exports: Tier-V (vertical column) and Tier-U (unstructured)
export half_levels, layer_thickness, pressure_coefficients
export cell_neighbor_table, cell_valence
# DUO icosahedral family (GRIDS_API §1 row 7, §10 loader-backed)
export DuoGrid, DuoLoader, build_duo_grid
export cell_centers, neighbors, metric_eval
export n_cells, n_vertices, n_edges, to_esm, family

# MPAS unstructured Voronoi family (GRIDS_API §1 row 6, §10 loader-backed)
export MpasGrid, MpasLoader, MpasMeshData, build_mpas_grid, mpas_mesh_data
export check_mesh, cell_area, edge_length, cell_distance, cell_center_cart
export max_edges

# Exports: Cartesian grid family (GRIDS_API.md §2.3)
export CartesianGrid, cell_widths, cell_volume

# Exports: Vertical grid family (GRIDS_API.md §2.3; column topology)
export VerticalGrid

# Exports: Lat-lon grid family (GRIDS_API.md §2.3; rectilinear sphere surface)
export LatLonGrid, nlon, nlon_uniform, row_offset, lon_edges, lon_centers, cell_area

# Exports: Connectivity
export EdgeDirection, West, East, South, North
export PanelNeighbor, PANEL_CONNECTIVITY
export transform_indices, verify_connectivity

# Exports: Metrics
export gnomonic_to_lonlat, gnomonic_to_cart
export gnomonic_metric, compute_cell_area, compute_edge_length
export compute_coord_jacobian, compute_forward_jacobian, compute_second_coord_jacobian
export tangent_vectors_3d, compute_edge_rotation_matrix

# Exports: Staggering
export VarLocation, CellCenter, UEdge, VEdge, Corner
export grid_size, full_array_size, ghost_array_size

# Exports: Discrete space and ghost cells
export CubedSphereDiscreteSpace, allocate_variable
export fill_ghost_cells!, extend_with_ghosts
export fill_ghost_cells_vector!, extend_with_ghosts_vector
export ghost_fill_indices, ghost_fill_arrayop

# Exports: ArrayOp utilities
export const_wrap, get_idx_vars, make_arrayop, evaluate_arrayop

# Exports: FV operators
export fv_divergence
export fv_gradient_xi, fv_gradient_eta
export fv_laplacian
export ppm_reconstruction, ppm_reconstruction!
export flux_1d
export transport_2d
export vertical_remap, vertical_remap_tendency

# Exports: Numerical transport operators
export flux_1d_ppm!, transport_2d_linrood!

# Exports: ArrayOp-based PPM transport operators
export flux_1d_ppm_arrayop, flux_to_tendency_arrayop, advective_tendency_arrayop
export compute_courant_numbers, compute_courant_numbers_arrayop
export transport_2d_ppm_arrayop
export ppm_reconstruction_arrayop

# Exports: FV3 super-grid
export compute_super_grid, compute_super_grid!, compute_angle_at_point

# Exports: FV3 wind operators
export covariant_to_contravariant, contravariant_to_covariant
export dgrid_to_cgrid!, dgrid_to_cgrid, dgrid_to_cgrid_arrayop
export compute_flux_with_sinsg!, compute_flux_with_sinsg
export compute_flux_with_sinsg_xi_arrayop, compute_flux_with_sinsg_eta_arrayop

# Exports: FV3 vorticity and kinetic energy
export fv_vorticity!, fv_vorticity, fv_vorticity_arrayop
export fv_vorticity_cellmean!, fv_vorticity_cellmean, fv_vorticity_cellmean_arrayop
export fv_absolute_vorticity!, fv_absolute_vorticity, fv_absolute_vorticity_arrayop
export fv_kinetic_energy!, fv_kinetic_energy, fv_kinetic_energy_arrayop
export fv_kinetic_energy_cell!, fv_kinetic_energy_cell, fv_kinetic_energy_cell_arrayop

# Exports: FV3 two-sided PPM
export ppm_edge_value_twosided, ppm_edge_value_twosided_limited
export flux_1d_ppm_twosided!

# Exports: Initial-condition projection
export project_initial_condition

# Exports: Rule catalog
export RuleFile, load_rules
export eval_coeff

# Exports: Arakawa staggering runtime
export ArakawaGrid, ArakawaStagger
export ArakawaA, ArakawaB, ArakawaC, ArakawaD, ArakawaE
export ArakawaBaseGrid, CartesianBase
export cell_centers, u_face, v_face, corners, neighbors, metric_eval
export arakawa_variable_locations, arakawa_location_shape, arakawa_shape, variable_shape
export to_esm

# -----------------------------------------------------------------------------
# `EarthSciDiscretizations.grids.<family>` generator namespace (GRIDS_API.md §2.3).
#
# Each family is exposed as a keyword-only function; the arakawa family is the
# first to land here. Future families (cartesian, lat_lon, cubed_sphere, …)
# add siblings in this submodule as their phase beads complete.
# -----------------------------------------------------------------------------
module grids

    using ..EarthSciDiscretizations: ArakawaGrid, ArakawaStagger, ArakawaBaseGrid,
        ArakawaA, ArakawaB, ArakawaC, ArakawaD, ArakawaE
    import ..EarthSciDiscretizations: _cartesian, _vertical, _latlon, build_duo_grid, DuoGrid, DuoLoader,
        build_mpas_grid, MpasGrid, MpasLoader, MpasMeshData, mpas_mesh_data

    """
        EarthSciDiscretizations.grids.duo(; loader, R=6.371e6, dtype=Float64, ghosts=0) -> DuoGrid

    Loader-backed icosahedral triangular grid (Heikes et al. 2023). Only
    `builtin://icosahedral/<level>` loader paths are honored today; `.duo`
    mesh files land with the ESS file-format spec.
    """
    duo(; kwargs...) = build_duo_grid(; kwargs...)

    """
        EarthSciDiscretizations.grids.mpas(; mesh=nothing, loader=nothing,
                                             reader_fn=nothing,
                                             R=6.371e6, dtype=Float64, ghosts=0)
            -> MpasGrid

    Loader-backed unstructured Voronoi grid (MPAS). Either `mesh`
    (`MpasMeshData`) or `loader` + `reader_fn` is required per
    GRIDS_API.md §10; NetCDF I/O is not bundled with this package.
    """
    mpas(; kwargs...) = build_mpas_grid(; kwargs...)

    """
        EarthSciDiscretizations.grids.cartesian(; nx, ny=nothing, nz=nothing,
                                                  extent=nothing, edges=nothing,
                                                  dtype=Float64, ghosts=0)
            -> CartesianGrid{T,N}

    Construct a 1D / 2D / 3D Cartesian grid (uniform via `extent` + `nx`[/`ny`/`nz`]
    or non-uniform via per-axis `edges`). See `GRIDS_API.md` §2.3 for the
    cross-binding API contract.
    """
    const cartesian = _cartesian

    """
        EarthSciDiscretizations.grids.vertical(; coordinate, nz=nothing, levels=nothing,
                                                 ak=nothing, bk=nothing, p0=1.0e5,
                                                 transition=nothing,
                                                 dtype=Float64, ghosts=0)
            -> VerticalGrid{T}

    Construct a 1D vertical column. Accepts coordinate kinds `:sigma`,
    `:eta`, `:z`, `:theta`, `:hybrid_sigma_theta`, `:z_star` (or the
    equivalent string literals). See `GRIDS_API.md` §2.3 for the
    cross-binding API contract and `src/grids/vertical.jl` for the
    coordinate-specific option table.
    """
    const vertical = _vertical

    """
        EarthSciDiscretizations.grids.lat_lon(; variant=:regular, nlon=nothing, nlat=nothing,
                                                nlon_per_row=nothing, lat_edges=nothing,
                                                lat_centers=nothing, R=6.371e6,
                                                dtype=Float64, ghosts=0,
                                                pole_policy=:none, lon_start=nothing)
            -> LatLonGrid{T}

    Construct a regular or reduced-Gaussian lat-lon grid on a sphere of radius
    `R`. Only `pole_policy=:none` is implemented in this phase. See
    `GRIDS_API.md` §2.3 and `src/grids/latlon.jl` for the option table and the
    cross-binding signature contract.
    """
    const lat_lon = _latlon

    """
        EarthSciDiscretizations.grids.arakawa(; base, stagger, ghosts=0, dtype=Float64)
            -> ArakawaGrid

    Construct an Arakawa-staggered grid over an underlying `base` grid. `stagger`
    accepts an `ArakawaStagger` value or the symbol literal `:A`, `:B`, `:C`,
    `:D`, `:E`. See `GRIDS_API.md` §2.3.
    """
    function arakawa(;
            base::Union{ArakawaBaseGrid, Nothing} = nothing,
            stagger::Union{ArakawaStagger, Symbol, Nothing} = nothing,
            ghosts::Int = 0,
            dtype::Type = Float64
        )
        base === nothing && throw(ArgumentError("arakawa: keyword argument `base` is required"))
        stagger === nothing && throw(ArgumentError("arakawa: keyword argument `stagger` is required"))
        s = stagger isa ArakawaStagger ? stagger : _parse_stagger(stagger)
        return ArakawaGrid(base, s; ghosts = ghosts, dtype = dtype)
    end

    function _parse_stagger(s::Symbol)
        return s === :A ? ArakawaA :
            s === :B ? ArakawaB :
            s === :C ? ArakawaC :
            s === :D ? ArakawaD :
            s === :E ? ArakawaE :
            throw(DomainError(s, "Unknown Arakawa stagger; expected :A, :B, :C, :D, or :E"))
    end

end # module grids

end # module
