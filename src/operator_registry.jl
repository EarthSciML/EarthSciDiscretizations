"""
Operator registry: maps operator type symbols to FV operators.
"""

struct OperatorRegistry
    grid::CubedSphereGrid
    disc::FVCubedSphere
end

build_operator_registry(grid::CubedSphereGrid, disc::FVCubedSphere) = OperatorRegistry(grid, disc)

function apply_operator(registry::OperatorRegistry, op_type::Symbol, args...)
    grid = registry.grid
    op_type == :gradient_xi && return fv_gradient_xi(args[1], grid)
    op_type == :gradient_eta && return fv_gradient_eta(args[1], grid)
    op_type == :divergence && return fv_divergence(args[1], args[2], grid)
    op_type == :laplacian && return fv_laplacian(args[1], grid)
    op_type == :flux_xi && return flux_1d(args[1], args[2], grid, :xi)
    op_type == :flux_eta && return flux_1d(args[1], args[2], grid, :eta)
    op_type == :transport_2d && return transport_2d(args[1], args[2], args[3], grid)
    error("Unknown operator type: $op_type")
end
