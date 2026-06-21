module EarthSciDiscretizations

using LinearAlgebra: norm, cross, dot
using SymbolicUtils: SymReal, BSImpl, idxs_for_arrayop, @arrayop, @syms
using SymbolicUtils
using Symbolics: unwrap, wrap, Differential, iscall, operation, arguments, Num
import Symbolics

# Grid infrastructure
include("grids/abstract_grid.jl")
include("grids/duo.jl")
include("grids/mpas.jl")
include("grids/cartesian.jl")
include("grids/vertical.jl")
include("grids/latlon.jl")

# Staggering and discrete space
include("staggering.jl")

# Arakawa staggering runtime (depends on VarLocation from staggering.jl)
include("grids/arakawa.jl")

# Discretization rule catalog (parser delegator; ESS integration pending)
include("rules.jl")

# AST coefficient evaluator (thin passthrough to EarthSciSerialization)
include("rule_eval.jl")

# DUO primal topology as value-invention FAQ — thin bridge routing edge
# enumeration / dedup / rank / inversion-join through the landed ESS M3 engine
# (EarthSciSerialization.Relational). Declarative companion:
# discretizations/grids/duo/primal_topology.esm.
include("topology_faq.jl")

# DUO grid CONSTRUCTION as the front-door declarative path (esd-ohd / W1):
# subdivision (D3), primal geometry (D2a), and edge/circumcenter geometry (D2b)
# bridges. `build_duo_grid` (grids/duo.jl) calls these + `primal_topology_faq`
# (D1a) instead of the imperative construction, routing every coordinate through
# `eval_coeff` and connectivity through `EarthSciSerialization.Relational`.
# Declarative companions under discretizations/grids/duo/{faq,rules}/*.esm.
include("subdivide_faq.jl")
include("primal_geometry_faq.jl")
include("edge_dual_geometry_faq.jl")

# MPAS-SCVT spherical-topology LEAF (esd-e5m.3 / D3): the one-time post-convergence
# topology op of the declarative SCVT mesh generator. The irreducible spherical
# Delaunay (convex hull of the generator directions) — canonical executor the
# s2bindings.rs S2B FFI (s2b-s7b) — composed with the landed declarative dual
# topology FAQ (voronoi_dual_topology_faq) to emit deterministic Voronoi
# connectivity. Determinism (integer) + tolerance (geometry) contract:
# discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md. Depends on
# topology_faq.jl (above).
include("grids/mpas_scvt_topology.jl")

# Cartesian grid construction as the structured-grid FAQ template (esd-3we.1 / S1)
# — thin bridge routing the affine coordinate map + elementwise metric/neighbor/
# boundary derivations through the landed ESS M1 evaluator (eval_coeff). Declarative
# companion: discretizations/grids/cartesian/rules/cartesian_construction.esm.
include("cartesian_faq.jl")

# Vertical (1-D column) grid construction as elementwise FAQ (esd-3we.2 / S2) —
# reuses the structured-grid template: a thin bridge routing the per-coordinate
# level synthesis + center/width/volume + named layer metrics through the landed
# ESS M1 evaluator (eval_coeff), with structural neighbor/boundary arrays. Declarative
# companion: discretizations/grids/vertical/rules/vertical_construction.esm.
include("vertical_faq.jl")

# Arakawa staggered-grid construction as the structured-grid FAQ (esd-3we.4 / S4) —
# thin bridge routing dx/dy + the two staggered 1-D affine maps + the cell-volume
# product through the landed ESS M1 evaluator; the four staggered locations, neighbor/
# boundary maps and static A–E location table are structural. Declarative companion:
# discretizations/grids/arakawa/rules/arakawa_construction.esm. Depends on the arakawa
# runtime (grids/arakawa.jl) and eval_coeff (rule_eval.jl), both included above.
include("arakawa_faq.jl")

# Lat-lon grid construction as the structured-grid FAQ (esd-3we.3 / S2) — thin bridge
# routing the per-row affine longitude map, the spherical-rectangle cell area, and the
# closed-form curvilinear (sin/cos) metric through the landed ESS M1 evaluator
# (eval_coeff); periodic / reduced-Gaussian nearest-center neighbor maps and boundary
# masks are structural index logic. Declarative companion:
# discretizations/grids/latlon/rules/latlon_construction.esm.
include("latlon_faq.jl")

# Top-level ESM → ODEProblem constructor (esd-3ck). Loads a PDE component
# .esm file and an optional GDD, runs the canonical ESS pipeline, and
# returns a SciMLBase.ODEProblem ready for the caller to solve.
include("ode_problem.jl")

# Exports: Grid types
export AbstractGrid
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

# Exports: Staggering
export VarLocation, CellCenter, UEdge, VEdge, Corner

# Exports: Rule catalog
export RuleFile, load_rules
export eval_coeff

# Exports: DUO primal topology value-invention FAQ bridge (esd-heg.1 / D1a)
export primal_topology_faq

# Exports: DUO→MPAS Voronoi-dual topology value-invention FAQ bridge (esd-heg.2 / D1b)
export voronoi_dual_topology_faq

# Exports: MPAS-SCVT spherical-topology leaf (esd-e5m.3 / D3) — the one-time
# post-convergence spherical Delaunay (wrap S2B) + deterministic Voronoi connectivity.
export scvt_spherical_delaunay, scvt_voronoi_connectivity

# Exports: DUO construction FAQ bridges wired into build_duo_grid (esd-ohd / W1)
export duo_subdivide_faq                              # D3 subdivision (esd-heg.3)
export primal_geometry_faq                            # D2a primal geometry (esd-heg.6)
export duo_circumcenter_geo_faq, duo_edge_length_faq, duo_cell_distance_faq  # D2b (esd-heg.7)

# Exports: Cartesian construction FAQ bridge — structured-grid template (esd-3we.1 / S1)
export cartesian_construction_faq

# Exports: Vertical construction FAQ bridge — reuses the template (esd-3we.2 / S2)
export vertical_construction_faq

# Exports: Arakawa staggered construction FAQ bridge (esd-3we.4 / S4)
export arakawa_construction_faq

# Exports: Lat-lon construction FAQ bridge — structured-grid template (esd-3we.3 / S2)
export latlon_construction_faq

# Exports: ESM → ODEProblem constructor (esd-3ck)
export build_ode_problem

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
# first to land here. Future families add siblings in this submodule as their
# phase beads complete.
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
