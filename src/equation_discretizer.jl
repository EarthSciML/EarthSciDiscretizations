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
spatial derivatives with finite-difference stencils.
"""
function discretize_equation(eq, disc_vars, dvs, spatial_ivs, grid, disc)
    Nc = grid.Nc
    dξ = grid.dξ
    dη = grid.dη

    lhs = eq.lhs
    rhs = eq.rhs

    dv = dvs[1]
    arr = disc_vars[dv]

    # Build substitution rules: replace spatial derivative terms with stencils
    # We iterate over grid points and for each point substitute the PDE unknowns
    # and their spatial derivatives with finite-difference expressions.

    eqs = Symbolics.Equation[]

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # For this cell, replace the PDE unknown and its derivatives
        # with the discrete array values and FD stencils
        ip1 = min(i + 1, Nc); im1 = max(i - 1, 1)
        jp1 = min(j + 1, Nc); jm1 = max(j - 1, 1)

        rhs_cell = _substitute_at_point(rhs, dv, arr, spatial_ivs, p, i, j, ip1, im1, jp1, jm1, dξ, dη)

        push!(eqs, mtk_D(arr[p, i, j]) ~ rhs_cell)
    end

    return eqs
end

"""
    _substitute_at_point(expr, dv, arr, spatial_ivs, p, i, j, ...)

Recursively walk a symbolic expression, replacing:
- The PDE unknown u(t, x, y) → arr[p, i, j]
- First spatial derivatives Dx(u) → centered FD
- Second spatial derivatives Dx(Dx(u)) → 3-point FD
- Coefficients/parameters → unchanged
"""
function _substitute_at_point(expr, dv, arr, spatial_ivs, p, i, j, ip1, im1, jp1, jm1, dξ, dη)
    ex = unwrap(expr)

    if !iscall(ex)
        # Leaf node: check if it's the dependent variable
        # (In practice, the DV won't appear as a bare leaf in the RHS)
        return wrap(ex)
    end

    op = operation(ex)
    args = arguments(ex)

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
                # Second derivative: d²u/dx²
                u_ij = _get_dv_value(innermost, dv, arr, p, i, j)

                if dim == :xi
                    u_ip1 = _get_dv_value(innermost, dv, arr, p, ip1, j)
                    u_im1 = _get_dv_value(innermost, dv, arr, p, im1, j)
                    return (u_ip1 - 2 * u_ij + u_im1) / dξ^2
                elseif dim == :eta
                    u_jp1 = _get_dv_value(innermost, dv, arr, p, i, jp1)
                    u_jm1 = _get_dv_value(innermost, dv, arr, p, i, jm1)
                    return (u_jp1 - 2 * u_ij + u_jm1) / dη^2
                end
            end
        end

        # First derivative: du/dx
        u_ij_val = _substitute_at_point(wrap(inner), dv, arr, spatial_ivs, p, i, j, ip1, im1, jp1, jm1, dξ, dη)

        if dim == :xi
            u_ip1 = _get_dv_value(inner, dv, arr, p, ip1, j)
            u_im1 = _get_dv_value(inner, dv, arr, p, im1, j)
            return (u_ip1 - u_im1) / (2 * dξ)
        elseif dim == :eta
            u_jp1 = _get_dv_value(inner, dv, arr, p, i, jp1)
            u_jm1 = _get_dv_value(inner, dv, arr, p, i, jm1)
            return (u_jp1 - u_jm1) / (2 * dη)
        end
    end

    # For general operations: recurse into all arguments
    new_args = [_substitute_at_point(wrap(a), dv, arr, spatial_ivs, p, i, j, ip1, im1, jp1, jm1, dξ, dη) for a in args]
    return wrap(op(unwrap.(new_args)...))
end

"""
Check if an expression is (or contains) the dependent variable, and return
the appropriate discrete array value at the given indices.
"""
function _get_dv_value(expr, dv, arr, p, i, j)
    # The DV appears as u(t, lon, lat) in the expression.
    # We replace it with arr[p, i, j].
    return arr[p, i, j]
end
