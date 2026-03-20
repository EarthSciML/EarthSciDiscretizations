module EarthSciDiscretizations

using LinearAlgebra: norm, cross, dot
using SymbolicUtils: SymReal, BSImpl, idxs_for_arrayop, @arrayop, @syms
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
