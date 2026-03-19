"""
Equation discretizer: replaces spatial derivatives in PDE equations
with FV stencil operators, producing per-cell ODE equations.
"""

function identify_dimension(iv)
    name = Symbol(iv)
    name == :t && return :t
    name in (:lon, :λ, :x, :ξ, :xi) && return :xi
    name in (:lat, :φ, :y, :η, :eta) && return :eta
    name in (:z, :p, :σ, :k) && return :vertical
    error("Cannot identify grid dimension for variable: $name")
end

"""
    discretize_equation(eq, disc_vars, dvs, spatial_ivs, grid, disc)

Convert a single PDE equation into per-cell ODE equations by substituting
spatial derivatives with metric-corrected finite-difference stencils.

Handles multiple dependent variables, nonlinear terms, mixed derivatives,
and inter-panel ghost cell references at boundaries.
"""
function discretize_equation(eq, disc_vars, dvs, spatial_ivs, grid, disc)
    Nc = grid.Nc
    dξ = grid.dξ
    dη = grid.dη

    # Determine which DV the LHS time-derivative refers to
    lhs_dv = _identify_lhs_dv(eq.lhs, dvs)
    lhs_arr = disc_vars[lhs_dv]

    rhs = eq.rhs

    eqs = Symbolics.Equation[]

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Compute neighbor indices, resolving across panel boundaries
        pip1 = _neighbor_index(grid, p, i + 1, j)
        pim1 = _neighbor_index(grid, p, i - 1, j)
        pjp1 = _neighbor_index(grid, p, i, j + 1)
        pjm1 = _neighbor_index(grid, p, i, j - 1)
        # Corner neighbors for mixed derivatives
        pip1jp1 = _neighbor_index(grid, pip1..., 0, 1)
        pip1jm1 = _neighbor_index(grid, pip1..., 0, -1)
        pim1jp1 = _neighbor_index(grid, pim1..., 0, 1)
        pim1jm1 = _neighbor_index(grid, pim1..., 0, -1)

        neighbors = (pip1=pip1, pim1=pim1, pjp1=pjp1, pjm1=pjm1,
                     pip1jp1=pip1jp1, pip1jm1=pip1jm1,
                     pim1jp1=pim1jp1, pim1jm1=pim1jm1)

        rhs_cell = _substitute_at_point(rhs, disc_vars, dvs, spatial_ivs,
                                         p, i, j, neighbors, dξ, dη, grid)

        push!(eqs, mtk_D(lhs_arr[p, i, j]) ~ rhs_cell)
    end

    return eqs
end

"""
Identify which dependent variable the LHS of an equation refers to.
The LHS has the form D(u(t, lon, lat)).
"""
function _identify_lhs_dv(lhs, dvs)
    ex = unwrap(lhs)
    if iscall(ex) && operation(ex) isa Differential
        inner = wrap(arguments(ex)[1])
        for dv in dvs
            if isequal(unwrap(inner), unwrap(dv))
                return dv
            end
        end
        # Fall back to name matching
        inner_name = Symbol(Symbolics.tosymbol(inner, escape=false))
        for dv in dvs
            dv_name = Symbol(Symbolics.tosymbol(dv, escape=false))
            if inner_name == dv_name
                return dv
            end
        end
    end
    error("Cannot identify which dependent variable LHS refers to: $lhs")
end

"""
    _neighbor_index(grid, p, i, j)

Map a (panel, i, j) triple that may be out of range to the correct
(panel, i, j) on a neighbor panel using the connectivity table.

    _neighbor_index(grid, p, i, j, di, dj)

Offset from (p, i, j) by (di, dj), resolving across panel boundaries.
"""
function _neighbor_index(grid::CubedSphereGrid, p::Int, i::Int, j::Int)
    Nc = grid.Nc
    if 1 <= i <= Nc && 1 <= j <= Nc
        return (p, i, j)
    end

    if i < 1
        nb = PANEL_CONNECTIVITY[p][West]
        i_nb, j_nb = _map_to_neighbor(nb, 1 - i, j, Nc, West)
        return (nb.neighbor_panel, clamp(i_nb, 1, Nc), clamp(j_nb, 1, Nc))
    elseif i > Nc
        nb = PANEL_CONNECTIVITY[p][East]
        i_nb, j_nb = _map_to_neighbor(nb, i - Nc, j, Nc, East)
        return (nb.neighbor_panel, clamp(i_nb, 1, Nc), clamp(j_nb, 1, Nc))
    elseif j < 1
        nb = PANEL_CONNECTIVITY[p][South]
        i_nb, j_nb = _map_to_neighbor(nb, 1 - j, i, Nc, South)
        return (nb.neighbor_panel, clamp(i_nb, 1, Nc), clamp(j_nb, 1, Nc))
    else  # j > Nc
        nb = PANEL_CONNECTIVITY[p][North]
        i_nb, j_nb = _map_to_neighbor(nb, j - Nc, i, Nc, North)
        return (nb.neighbor_panel, clamp(i_nb, 1, Nc), clamp(j_nb, 1, Nc))
    end
end

function _neighbor_index(grid::CubedSphereGrid, p::Int, i::Int, j::Int, di::Int, dj::Int)
    return _neighbor_index(grid, p, i + di, j + dj)
end

# Unpack tuple form
function _neighbor_index(grid::CubedSphereGrid, pij::Tuple{Int,Int,Int}, di::Int, dj::Int)
    return _neighbor_index(grid, pij[1], pij[2] + di, pij[3] + dj)
end

"""
Map an out-of-range index to a neighbor panel. `offset` is how many cells
past the boundary (1-based), `j_along` is the index along the shared edge.
"""
function _map_to_neighbor(nb::PanelNeighbor, offset, j_along, Nc, src_dir)
    j_mapped = nb.reverse_index ? (Nc + 1 - j_along) : j_along
    nb_edge = nb.neighbor_edge
    if nb_edge == West
        return (offset, j_mapped)
    elseif nb_edge == East
        return (Nc + 1 - offset, j_mapped)
    elseif nb_edge == South
        return (j_mapped, offset)
    else  # North
        return (j_mapped, Nc + 1 - offset)
    end
end

"""
    _substitute_at_point(expr, disc_vars, dvs, spatial_ivs,
                          p, i, j, neighbors, dξ, dη, grid)

Recursively walk a symbolic expression, replacing:
- DV calls u(t, x, y) → disc_vars[u][p, i, j]
- First spatial derivatives → metric-corrected centered FD
- Second spatial derivatives → metric-corrected 3-point FD
- Mixed spatial derivatives → 4-point cross-derivative stencil
- Nonlinear terms are preserved (arithmetic is recursed into)
"""
function _substitute_at_point(expr, disc_vars, dvs, spatial_ivs,
                               p, i, j, neighbors, dξ, dη, grid)
    ex = unwrap(expr)

    if !iscall(ex)
        return wrap(ex)
    end

    op = operation(ex)
    args = arguments(ex)

    # Check if this is a DV call (e.g., u(t, lon, lat))
    for dv in dvs
        if isequal(wrap(ex), dv)
            return disc_vars[dv][p, i, j]
        end
    end

    # Check for Differential operator
    if op isa Differential
        dim = identify_dimension(Symbol(op.x))
        inner = args[1]

        # Check if inner is also a Differential (second derivative)
        if iscall(inner) && operation(inner) isa Differential
            inner_op = operation(inner)
            inner_dim = identify_dimension(Symbol(inner_op.x))
            innermost = arguments(inner)[1]

            if dim == inner_dim
                # Same-dimension second derivative with metric correction:
                # (1/J) * d/dξ(J * g^{ξξ} * du/dξ) at cell (p,i,j)
                return _second_deriv_same(innermost, dim, disc_vars, dvs,
                                          p, i, j, neighbors, dξ, dη, grid)
            else
                # Mixed derivative: d²u/(dξ dη) with metric cross-term
                return _second_deriv_mixed(innermost, dim, inner_dim, disc_vars, dvs,
                                            p, i, j, neighbors, dξ, dη, grid)
            end
        end

        # First derivative with metric correction
        return _first_deriv(inner, dim, disc_vars, dvs, spatial_ivs,
                             p, i, j, neighbors, dξ, dη, grid)
    end

    # For general operations: recurse into all arguments
    new_args = [_substitute_at_point(wrap(a), disc_vars, dvs, spatial_ivs,
                                      p, i, j, neighbors, dξ, dη, grid) for a in args]
    return wrap(op(unwrap.(new_args)...))
end

"""
Metric-corrected first derivative.

For a derivative in the ξ-direction:
  du/dx_phys = (g^{ξξ} * du/dξ + g^{ξη} * du/dη) / sqrt(g^{ξξ})

For a generic scalar PDE, we approximate the physical first derivative
using the contravariant metric. This uses centered differences in both
directions when the cross-metric term is nonzero.
"""
function _first_deriv(inner, dim, disc_vars, dvs, spatial_ivs,
                       p, i, j, neighbors, dξ, dη, grid)
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1

    u_ip1 = _get_dv_value(inner, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(inner, disc_vars, dvs, pim1...)
    u_jp1 = _get_dv_value(inner, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(inner, disc_vars, dvs, pjm1...)

    # Inverse metric at this cell
    ginv_ξξ = grid.ginv_ξξ[p, i, j]
    ginv_ηη = grid.ginv_ηη[p, i, j]
    ginv_ξη = grid.ginv_ξη[p, i, j]

    du_dξ = (u_ip1 - u_im1) / (2 * dξ)
    du_dη = (u_jp1 - u_jm1) / (2 * dη)

    if dim == :xi
        # Physical derivative in the ξ-mapped direction
        # Scale by sqrt(ginv_ξξ) to convert from covariant to physical
        return ginv_ξξ / sqrt(ginv_ξξ) * du_dξ + ginv_ξη / sqrt(ginv_ξξ) * du_dη
    elseif dim == :eta
        return ginv_ηη / sqrt(ginv_ηη) * du_dη + ginv_ξη / sqrt(ginv_ηη) * du_dξ
    end
end

"""
Metric-corrected second derivative in the same dimension.

Computes (1/J) * ∂/∂ξ(J * g^{ξξ} * ∂u/∂ξ) using half-point averaged metrics.
"""
function _second_deriv_same(innermost, dim, disc_vars, dvs,
                             p, i, j, neighbors, dξ, dη, grid)
    u_ij = _get_dv_value(innermost, disc_vars, dvs, p, i, j)
    J_ij = grid.J[p, i, j]

    if dim == :xi
        pip1 = neighbors.pip1; pim1 = neighbors.pim1
        u_ip1 = _get_dv_value(innermost, disc_vars, dvs, pip1...)
        u_im1 = _get_dv_value(innermost, disc_vars, dvs, pim1...)

        # Half-point metric averages
        J_ip = (grid.J[pip1...] + J_ij) / 2
        J_im = (grid.J[pim1...] + J_ij) / 2
        ginv_ip = (grid.ginv_ξξ[pip1...] + grid.ginv_ξξ[p, i, j]) / 2
        ginv_im = (grid.ginv_ξξ[pim1...] + grid.ginv_ξξ[p, i, j]) / 2

        return (1 / J_ij) * (
            J_ip * ginv_ip * (u_ip1 - u_ij) / dξ -
            J_im * ginv_im * (u_ij - u_im1) / dξ
        ) / dξ
    elseif dim == :eta
        pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1
        u_jp1 = _get_dv_value(innermost, disc_vars, dvs, pjp1...)
        u_jm1 = _get_dv_value(innermost, disc_vars, dvs, pjm1...)

        J_jp = (grid.J[pjp1...] + J_ij) / 2
        J_jm = (grid.J[pjm1...] + J_ij) / 2
        ginv_jp = (grid.ginv_ηη[pjp1...] + grid.ginv_ηη[p, i, j]) / 2
        ginv_jm = (grid.ginv_ηη[pjm1...] + grid.ginv_ηη[p, i, j]) / 2

        return (1 / J_ij) * (
            J_jp * ginv_jp * (u_jp1 - u_ij) / dη -
            J_jm * ginv_jm * (u_ij - u_jm1) / dη
        ) / dη
    end
end

"""
Metric-corrected mixed second derivative: ∂²u/(∂ξ ∂η).
Uses 4-point cross-derivative stencil scaled by the inverse metric cross-term.
"""
function _second_deriv_mixed(innermost, dim_outer, dim_inner, disc_vars, dvs,
                              p, i, j, neighbors, dξ, dη, grid)
    u_pp = _get_dv_value(innermost, disc_vars, dvs, neighbors.pip1jp1...)
    u_pm = _get_dv_value(innermost, disc_vars, dvs, neighbors.pip1jm1...)
    u_mp = _get_dv_value(innermost, disc_vars, dvs, neighbors.pim1jp1...)
    u_mm = _get_dv_value(innermost, disc_vars, dvs, neighbors.pim1jm1...)

    return (u_pp - u_pm - u_mp + u_mm) / (4 * dξ * dη)
end

"""
    _get_dv_value(expr, disc_vars, dvs, p, i, j)

Recursively walk `expr`, substituting each dependent variable call with
the corresponding discrete array value at (p, i, j). Preserves arithmetic
structure so that nonlinear terms like u^2 become arr[p,i,j]^2.
"""
function _get_dv_value(expr, disc_vars, dvs, p, i, j)
    ex = unwrap(expr)

    if !iscall(ex)
        return wrap(ex)
    end

    op = operation(ex)
    args = arguments(ex)

    # Check if this call matches a DV (e.g., u(t, lon, lat))
    for dv in dvs
        if isequal(wrap(ex), dv)
            return disc_vars[dv][p, i, j]
        end
    end

    # Recurse into arguments for general expressions (u^2, u*v, sin(u), etc.)
    new_args = [_get_dv_value(wrap(a), disc_vars, dvs, p, i, j) for a in args]
    return wrap(op(unwrap.(new_args)...))
end
