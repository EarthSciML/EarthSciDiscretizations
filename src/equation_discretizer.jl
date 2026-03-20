"""
Equation discretizer: replaces spatial derivatives in PDE equations
with FV stencil operators, producing per-cell ODE equations.
"""

function identify_dimension(iv)
    name = Symbol(iv)
    name == :t && return :t
    name in (:lon, :λ) && return :lon
    name in (:lat, :φ) && return :lat
    name in (:x, :ξ, :xi) && return :xi
    name in (:y, :η, :eta) && return :eta
    name in (:z, :p, :σ, :k) && return :vertical
    error("Cannot identify grid dimension for variable: $name")
end

"""
Return true if `dim` is a physical (lon/lat) coordinate that needs
chain-rule transformation to computational (ξ/η) coordinates.
"""
_is_physical_coord(dim::Symbol) = dim in (:lon, :lat)

"""
    discretize_equation(eq, disc_vars, dvs, spatial_ivs, grid, disc)

Convert a single PDE equation into per-cell ODE equations by substituting
spatial derivatives with metric-corrected finite-difference stencils.

Handles multiple dependent variables, nonlinear terms, mixed derivatives,
and inter-panel ghost cell references at boundaries.

Physical coordinates (lon, lat) are transformed to computational coordinates
(ξ, η) via the chain rule using the precomputed coordinate Jacobian.
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
- First spatial derivatives → chain-rule-transformed centered FD
- Second spatial derivatives → chain-rule-transformed 3-point FD
- Mixed spatial derivatives → chain-rule-transformed cross-derivative
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

            return _second_deriv(innermost, dim, inner_dim, disc_vars, dvs,
                                 p, i, j, neighbors, dξ, dη, grid)
        end

        # First derivative
        return _first_deriv(inner, dim, disc_vars, dvs, spatial_ivs,
                             p, i, j, neighbors, dξ, dη, grid)
    end

    # For general operations: recurse into all arguments
    new_args = [_substitute_at_point(wrap(a), disc_vars, dvs, spatial_ivs,
                                      p, i, j, neighbors, dξ, dη, grid) for a in args]
    return wrap(op(unwrap.(new_args)...))
end

"""
First derivative with chain-rule coordinate transformation.

For physical coordinates (lon, lat), applies:
    ∂u/∂lon = (∂ξ/∂lon)·∂u/∂ξ + (∂η/∂lon)·∂u/∂η
    ∂u/∂lat = (∂ξ/∂lat)·∂u/∂ξ + (∂η/∂lat)·∂u/∂η

For computational coordinates (ξ, η), uses direct centered differences.
"""
function _first_deriv(inner, dim, disc_vars, dvs, spatial_ivs,
                       p, i, j, neighbors, dξ, dη, grid)
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1

    u_ip1 = _get_dv_value(inner, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(inner, disc_vars, dvs, pim1...)
    u_jp1 = _get_dv_value(inner, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(inner, disc_vars, dvs, pjm1...)

    du_dξ = (u_ip1 - u_im1) / (2 * dξ)
    du_dη = (u_jp1 - u_jm1) / (2 * dη)

    if dim == :lon
        # Chain rule: d/dlon = (dξ/dlon)·d/dξ + (dη/dlon)·d/dη
        dξ_dlon = grid.dξ_dlon[p, i, j]
        dη_dlon = grid.dη_dlon[p, i, j]
        return dξ_dlon * du_dξ + dη_dlon * du_dη
    elseif dim == :lat
        # Chain rule: d/dlat = (dξ/dlat)·d/dξ + (dη/dlat)·d/dη
        dξ_dlat = grid.dξ_dlat[p, i, j]
        dη_dlat = grid.dη_dlat[p, i, j]
        return dξ_dlat * du_dξ + dη_dlat * du_dη
    elseif dim == :xi
        return du_dξ
    elseif dim == :eta
        return du_dη
    end
end

"""
Second derivative with full chain-rule coordinate transformation.

For physical coordinates, applies the complete chain rule including the
second-derivative correction terms that arise from the nonlinear mapping:
    ∂²u/∂lon² = (∂ξ/∂lon)²·∂²u/∂ξ² + 2(∂ξ/∂lon)(∂η/∂lon)·∂²u/(∂ξ∂η) + (∂η/∂lon)²·∂²u/∂η²
              + (∂²ξ/∂lon²)·∂u/∂ξ + (∂²η/∂lon²)·∂u/∂η

Uses plain centered second differences (not metric-weighted operators)
for the ∂²u/∂ξ² and ∂²u/∂η² terms, as required by the chain rule.

For computational coordinates, uses direct stencils.
"""
function _second_deriv(innermost, dim_outer, dim_inner, disc_vars, dvs,
                        p, i, j, neighbors, dξ, dη, grid)
    # Plain (non-metric-weighted) second differences in computational coords
    d2u_dξ2 = _plain_second_deriv_xi(innermost, disc_vars, dvs, p, i, j, neighbors, dξ)
    d2u_dη2 = _plain_second_deriv_eta(innermost, disc_vars, dvs, p, i, j, neighbors, dη)
    d2u_dξdη = _comp_mixed_deriv(innermost, disc_vars, dvs, p, i, j, neighbors, dξ, dη)

    # First derivatives (needed for second-order correction terms)
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1
    u_ip1 = _get_dv_value(innermost, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(innermost, disc_vars, dvs, pim1...)
    u_jp1 = _get_dv_value(innermost, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(innermost, disc_vars, dvs, pjm1...)
    du_dξ = (u_ip1 - u_im1) / (2 * dξ)
    du_dη = (u_jp1 - u_jm1) / (2 * dη)

    # Chain-rule coefficients (first derivatives of coordinate map)
    a_ξ_outer, a_η_outer = _chain_rule_coeffs(dim_outer, grid, p, i, j)
    a_ξ_inner, a_η_inner = _chain_rule_coeffs(dim_inner, grid, p, i, j)

    # Second-derivative chain-rule correction coefficients: (∂²ξ/∂x∂y, ∂²η/∂x∂y)
    b_ξ, b_η = _second_chain_rule_coeffs(dim_outer, dim_inner, grid, p, i, j)

    # Full chain rule for second derivatives:
    # ∂²u/∂x∂y = Σ_k Σ_l (∂ξ_k/∂x)(∂ξ_l/∂y) ∂²u/(∂ξ_k ∂ξ_l)
    #           + Σ_k (∂²ξ_k/∂x∂y) ∂u/∂ξ_k
    return (a_ξ_outer * a_ξ_inner * d2u_dξ2 +
            (a_ξ_outer * a_η_inner + a_η_outer * a_ξ_inner) * d2u_dξdη +
            a_η_outer * a_η_inner * d2u_dη2 +
            b_ξ * du_dξ + b_η * du_dη)
end

"""
Return (coeff_ξ, coeff_η) for the chain rule transformation of dimension `dim`.
For computational coordinates: (:xi → (1,0), :eta → (0,1)).
For physical coordinates: uses the coordinate Jacobian.
"""
function _chain_rule_coeffs(dim, grid, p, i, j)
    if dim == :lon
        return (grid.dξ_dlon[p, i, j], grid.dη_dlon[p, i, j])
    elseif dim == :lat
        return (grid.dξ_dlat[p, i, j], grid.dη_dlat[p, i, j])
    elseif dim == :xi
        return (1.0, 0.0)
    elseif dim == :eta
        return (0.0, 1.0)
    end
end

"""
Return (coeff_ξ, coeff_η) for the second-derivative chain-rule correction.
These are the (∂²ξ/∂x∂y, ∂²η/∂x∂y) terms that arise because the coordinate
mapping is nonlinear.
"""
function _second_chain_rule_coeffs(dim_outer, dim_inner, grid, p, i, j)
    if dim_outer == :lon && dim_inner == :lon
        return (grid.d2ξ_dlon2[p, i, j], grid.d2η_dlon2[p, i, j])
    elseif dim_outer == :lat && dim_inner == :lat
        return (grid.d2ξ_dlat2[p, i, j], grid.d2η_dlat2[p, i, j])
    elseif (dim_outer == :lon && dim_inner == :lat) || (dim_outer == :lat && dim_inner == :lon)
        return (grid.d2ξ_dlondlat[p, i, j], grid.d2η_dlondlat[p, i, j])
    else
        # Computational coordinates: no second-order correction
        return (0.0, 0.0)
    end
end

"""
Plain centered second difference in ξ: (u_{i+1} - 2u_i + u_{i-1}) / dξ².
No metric weighting — this is the form required by the chain rule.
"""
function _plain_second_deriv_xi(innermost, disc_vars, dvs,
                                 p, i, j, neighbors, dξ)
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    u_ij = _get_dv_value(innermost, disc_vars, dvs, p, i, j)
    u_ip1 = _get_dv_value(innermost, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(innermost, disc_vars, dvs, pim1...)
    return (u_ip1 - 2 * u_ij + u_im1) / (dξ^2)
end

"""
Plain centered second difference in η: (u_{j+1} - 2u_j + u_{j-1}) / dη².
No metric weighting — this is the form required by the chain rule.
"""
function _plain_second_deriv_eta(innermost, disc_vars, dvs,
                                  p, i, j, neighbors, dη)
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1
    u_ij = _get_dv_value(innermost, disc_vars, dvs, p, i, j)
    u_jp1 = _get_dv_value(innermost, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(innermost, disc_vars, dvs, pjm1...)
    return (u_jp1 - 2 * u_ij + u_jm1) / (dη^2)
end

"""
Covariant second derivative in ξ:
    (1/J)·∂/∂ξ(J·g^{ξξ}·∂u/∂ξ) with half-point averaged metrics.
Used by the covariant Laplacian (not the chain-rule path).
"""
function _covariant_second_deriv_xi(innermost, disc_vars, dvs,
                                     p, i, j, neighbors, dξ, grid)
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    u_ij = _get_dv_value(innermost, disc_vars, dvs, p, i, j)
    u_ip1 = _get_dv_value(innermost, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(innermost, disc_vars, dvs, pim1...)

    J_ij = grid.J[p, i, j]
    J_ip = (grid.J[pip1...] + J_ij) / 2
    J_im = (grid.J[pim1...] + J_ij) / 2
    ginv_ip = (grid.ginv_ξξ[pip1...] + grid.ginv_ξξ[p, i, j]) / 2
    ginv_im = (grid.ginv_ξξ[pim1...] + grid.ginv_ξξ[p, i, j]) / 2

    return (1 / J_ij) * (
        J_ip * ginv_ip * (u_ip1 - u_ij) / dξ -
        J_im * ginv_im * (u_ij - u_im1) / dξ
    ) / dξ
end

"""
Covariant second derivative in η:
    (1/J)·∂/∂η(J·g^{ηη}·∂u/∂η) with half-point averaged metrics.
Used by the covariant Laplacian (not the chain-rule path).
"""
function _covariant_second_deriv_eta(innermost, disc_vars, dvs,
                                      p, i, j, neighbors, dη, grid)
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1
    u_ij = _get_dv_value(innermost, disc_vars, dvs, p, i, j)
    u_jp1 = _get_dv_value(innermost, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(innermost, disc_vars, dvs, pjm1...)

    J_ij = grid.J[p, i, j]
    J_jp = (grid.J[pjp1...] + J_ij) / 2
    J_jm = (grid.J[pjm1...] + J_ij) / 2
    ginv_jp = (grid.ginv_ηη[pjp1...] + grid.ginv_ηη[p, i, j]) / 2
    ginv_jm = (grid.ginv_ηη[pjm1...] + grid.ginv_ηη[p, i, j]) / 2

    return (1 / J_ij) * (
        J_jp * ginv_jp * (u_jp1 - u_ij) / dη -
        J_jm * ginv_jm * (u_ij - u_jm1) / dη
    ) / dη
end

"""
Covariant cross-derivative term for the Laplacian:
    (1/J)·∂/∂ξ(J·g^{ξη}·∂u/∂η) + (1/J)·∂/∂η(J·g^{ξη}·∂u/∂ξ)

Expands to:
    2·g^{ξη}·∂²u/(∂ξ∂η) + (1/J)·[∂(J·g^{ξη})/∂ξ·∂u/∂η + ∂(J·g^{ξη})/∂η·∂u/∂ξ]
"""
function _covariant_cross_deriv(innermost, disc_vars, dvs,
                                 p, i, j, neighbors, dξ, dη, grid)
    d2u_dξdη = _comp_mixed_deriv(innermost, disc_vars, dvs, p, i, j, neighbors, dξ, dη)

    # First derivatives of u
    pip1 = neighbors.pip1; pim1 = neighbors.pim1
    pjp1 = neighbors.pjp1; pjm1 = neighbors.pjm1
    u_ip1 = _get_dv_value(innermost, disc_vars, dvs, pip1...)
    u_im1 = _get_dv_value(innermost, disc_vars, dvs, pim1...)
    u_jp1 = _get_dv_value(innermost, disc_vars, dvs, pjp1...)
    u_jm1 = _get_dv_value(innermost, disc_vars, dvs, pjm1...)
    du_dξ = (u_ip1 - u_im1) / (2 * dξ)
    du_dη = (u_jp1 - u_jm1) / (2 * dη)

    # Derivatives of J·g^{ξη} using centered differences
    Jgxe_ip1 = grid.J[pip1...] * grid.ginv_ξη[pip1...]
    Jgxe_im1 = grid.J[pim1...] * grid.ginv_ξη[pim1...]
    Jgxe_jp1 = grid.J[pjp1...] * grid.ginv_ξη[pjp1...]
    Jgxe_jm1 = grid.J[pjm1...] * grid.ginv_ξη[pjm1...]
    dJgxe_dξ = (Jgxe_ip1 - Jgxe_im1) / (2 * dξ)
    dJgxe_dη = (Jgxe_jp1 - Jgxe_jm1) / (2 * dη)

    J_center = grid.J[p, i, j]
    first_deriv_correction = (1 / J_center) * (dJgxe_dξ * du_dη + dJgxe_dη * du_dξ)

    return 2 * grid.ginv_ξη[p, i, j] * d2u_dξdη + first_deriv_correction
end

"""
Computational-coordinate mixed derivative ∂²u/(∂ξ∂η) using 4-point stencil.
"""
function _comp_mixed_deriv(innermost, disc_vars, dvs,
                            p, i, j, neighbors, dξ, dη)
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
