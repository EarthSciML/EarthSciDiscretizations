module EarthSciDiscretizations

using LinearAlgebra: norm, cross, dot
using SymbolicUtils: SymReal, BSImpl, idxs_for_arrayop, @arrayop, @makearray, @syms
using SymbolicUtils
using Symbolics: unwrap, wrap, Differential, iscall, operation, arguments, Num
import Symbolics
using ModelingToolkit: PDESystem, System, mtkcompile, ODEProblem, @named
using ModelingToolkit: t_nounits as mtk_t, D_nounits as mtk_D
import ModelingToolkit
import SciMLBase
import DomainSets

# Grid infrastructure
include("grids/abstract_grid.jl")
include("grids/panel_connectivity.jl")
include("grids/metric_tensors.jl")
include("grids/super_grid.jl")
include("grids/cubed_sphere.jl")

# Staggering and discrete space
include("staggering.jl")
include("discrete_space.jl")
include("ghost_cells.jl")

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

# Precomputed FV stencils
include("fv_stencil.jl")

# Discretization pipeline
include("discretization.jl")
include("equation_discretizer.jl")
include("bc_handler.jl")

# Exports: Grid types
export AbstractGrid, AbstractCubedSphereGrid, CubedSphereGrid
export total_area

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

# Exports: Precomputed stencils
export FVLaplacianStencil, FVGradientStencil
export precompute_laplacian_stencil, precompute_gradient_stencil
export apply_laplacian!, apply_gradient!

# Exports: Discretization pipeline
export FVCubedSphere
export fv_laplacian_extended, fv_gradient_extended
export project_initial_condition
export identify_dimension

end # module
