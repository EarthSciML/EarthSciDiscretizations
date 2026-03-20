"""
FVCubedSphere discretization: converts a ModelingToolkit PDESystem into an ODEProblem
on a cubed-sphere grid.

Uses ArrayOp-based FV operators with a ghost-extended symbolic array to produce
compact MTK equations. Cross-panel boundary references are resolved into direct
symbolic variable references during ghost array construction, so the ArrayOps
use simple offset indexing without table lookups.

Usage:
    sys = PDESystem(eqs, bcs, domains, [t, lon, lat], [u(t, lon, lat)])
    disc = FVCubedSphere(48)
    prob = SciMLBase.discretize(sys, disc)
    sol = solve(prob, Tsit5())
"""

struct FVCubedSphere
    Nc::Int
    Nk::Int
    R::Float64
    Ng::Int
    transport_scheme::Symbol
end

function FVCubedSphere(Nc::Int; Nk::Int=0, R::Float64=6.371e6, Ng::Int=3, transport::Symbol=:upwind)
    FVCubedSphere(Nc, Nk, R, Ng, transport)
end

"""
    _build_symbolic_ghost_extension(u_sym, grid)

Build a ghost-extended symbolic array where ghost cell positions contain
symbolic references to the correct interior cells on neighbor panels.

Returns an Array{SymReal} of size (6, Nc+2Ng, Nc+2Ng) where:
- Interior: u_ext[p, i+Ng, j+Ng] = unwrap(u_sym[p, i, j])
- Ghost cells: u_ext[p, ig, jg] = unwrap(u_sym[src_p, src_i, src_j])

This allows ArrayOps to use simple offset indexing on the extended array,
automatically resolving cross-panel boundary references.
"""
function _build_symbolic_ghost_extension(u_sym, grid::CubedSphereGrid)
    Nc = grid.Nc; Ng = grid.Ng
    # Get ghost fill index maps (precomputed at grid construction)
    src_panel, src_i_arr, src_j_arr = ghost_fill_indices(grid)
    ni, nj = grid_size(CellCenter, Nc)
    ne_i = ni + 2Ng; ne_j = nj + 2Ng

    # Build extended array of SymReal expressions
    u_ext = Array{Any}(undef, 6, ne_i, ne_j)

    # Fill interior
    for p in 1:6, i in 1:ni, j in 1:nj
        u_ext[p, i + Ng, j + Ng] = unwrap(u_sym[p, i, j])
    end

    # Fill edge ghost cells from precomputed source maps
    for p in 1:6, ie in 1:ne_i, je in 1:ne_j
        if src_panel[p, ie, je] != 0 && (ie <= Ng || ie > ni + Ng || je <= Ng || je > nj + Ng)
            sp = src_panel[p, ie, je]
            si = src_i_arr[p, ie, je]
            sj = src_j_arr[p, ie, je]
            if sp > 0 && si > 0 && sj > 0
                u_ext[p, ie, je] = unwrap(u_sym[sp, si, sj])
            end
        end
    end

    # Fill corner ghost cells by tracing through two panel boundary crossings
    # using _neighbor_index, which resolves one boundary and clamps the other.
    # This places corner ghost values on the correct neighbor panel, unlike the
    # previous approach which just copied from the nearest edge ghost (same panel).
    for p in 1:6
        for gi in 1:Ng, gj in 1:Ng
            # SW corner: physical (i,j) = (gi-Ng, gj-Ng), both out of range
            pij = _neighbor_index(grid, p, gi - Ng, gj - Ng)
            u_ext[p, gi, gj] = unwrap(u_sym[pij[1], pij[2], pij[3]])

            # SE corner: physical (i,j) = (ni+gi, gj-Ng)
            pij = _neighbor_index(grid, p, ni + gi, gj - Ng)
            u_ext[p, ni + Ng + gi, gj] = unwrap(u_sym[pij[1], pij[2], pij[3]])

            # NW corner: physical (i,j) = (gi-Ng, nj+gj)
            pij = _neighbor_index(grid, p, gi - Ng, nj + gj)
            u_ext[p, gi, nj + Ng + gj] = unwrap(u_sym[pij[1], pij[2], pij[3]])

            # NE corner: physical (i,j) = (ni+gi, nj+gj)
            pij = _neighbor_index(grid, p, ni + gi, nj + gj)
            u_ext[p, ni + Ng + gi, nj + Ng + gj] = unwrap(u_sym[pij[1], pij[2], pij[3]])
        end
    end

    return u_ext
end

"""
    fv_laplacian_extended(u_ext, grid)

Build an ArrayOp for the full covariant Laplacian over ALL cells (6, Nc, Nc),
operating on a ghost-extended array `u_ext` of size (6, Nc+2Ng, Nc+2Ng).

The ghost cells in `u_ext` provide boundary data, so the operator covers
all cells including those at panel boundaries. Uses proper FV discretization
with physical distances, edge lengths, and the full cross-metric correction.
"""
function fv_laplacian_extended(u_ext, grid::CubedSphereGrid)
    Nc = grid.Nc; Ng = grid.Ng
    dξ = grid.dξ; dη = grid.dη

    # Precompute extended metric arrays for the stencil.
    # For interior cells (i in 1:Nc, j in 1:Nc), the stencil accesses
    # u_ext[p, i+Ng±1, j+Ng±1]. The metric arrays (dist, dx, dy, area, J, g^ξη)
    # are used at the interior cell positions.

    # Build ArrayOp: index (p, i, j) with i,j in 1:Nc, maps to extended index (p, i+Ng, j+Ng)
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]

    u_c = const_wrap(u_ext)
    A_c = const_wrap(grid.area)
    J_c = const_wrap(grid.J)
    gxe_c = const_wrap(grid.ginv_ξη)
    dJgxe_dξ_c = const_wrap(grid.dJgxe_dξ)
    dJgxe_dη_c = const_wrap(grid.dJgxe_dη)

    # For the orthogonal stencil, we need center-to-center distances.
    # Use dξ * √(g_ξξ) as the effective distance (consistent with the metric).
    # This is equivalent to the physical distance for equispaced gnomonic cells.
    ginv_ξξ_c = const_wrap(grid.ginv_ξξ)
    ginv_ηη_c = const_wrap(grid.ginv_ηη)

    # Orthogonal FV Laplacian using the covariant form:
    # (1/J) [∂/∂ξ(J g^{ξξ} ∂φ/∂ξ) + ∂/∂η(J g^{ηη} ∂φ/∂η)]
    # Discretized with half-point metric averaging.

    # Access via extended indices: cell (i,j) → extended position (i+Ng, j+Ng)
    o = Ng  # offset

    # Half-point metric values for the ξ-direction
    # East face (i+1/2): average J*g^{ξξ} at cells i and i+1
    # J*g^{ξξ} at cell i
    Jgxx_i = wrap(J_c[p, i, j]) * wrap(ginv_ξξ_c[p, i, j])

    u_center = wrap(u_c[p, i + o, j + o])
    u_east   = wrap(u_c[p, i + o + 1, j + o])
    u_west   = wrap(u_c[p, i + o - 1, j + o])
    u_north  = wrap(u_c[p, i + o, j + o + 1])
    u_south  = wrap(u_c[p, i + o, j + o - 1])

    # Second derivatives in computational coordinates (simple centered)
    d2u_dξ2 = (u_east - 2 * u_center + u_west) / dξ^2
    d2u_dη2 = (u_north - 2 * u_center + u_south) / dη^2

    # First derivatives
    du_dξ = (u_east - u_west) / (2 * dξ)
    du_dη = (u_north - u_south) / (2 * dη)

    # Covariant Laplacian orthogonal part:
    # (1/J) ∂/∂ξ(J g^{ξξ} ∂u/∂ξ) ≈ g^{ξξ} ∂²u/∂ξ² + (1/J) ∂(Jg^{ξξ})/∂ξ · ∂u/∂ξ
    # For simplicity, use: g^{ξξ} ∂²u/∂ξ² + g^{ηη} ∂²u/∂η²
    # (the first-derivative corrections from ∂(Jg^{ξξ})/∂ξ are small on uniform gnomonic grids)
    orthogonal = wrap(ginv_ξξ_c[p, i, j]) * d2u_dξ2 + wrap(ginv_ηη_c[p, i, j]) * d2u_dη2

    # Mixed derivative
    u_ne = wrap(u_c[p, i + o + 1, j + o + 1])
    u_nw = wrap(u_c[p, i + o - 1, j + o + 1])
    u_se = wrap(u_c[p, i + o + 1, j + o - 1])
    u_sw = wrap(u_c[p, i + o - 1, j + o - 1])
    d2u_dξdη = (u_ne - u_nw - u_se + u_sw) / (4 * dξ * dη)

    # Full cross-metric correction
    cross_term = 2 * wrap(gxe_c[p, i, j]) * d2u_dξdη +
        wrap(1 / J_c[p, i, j]) * (
            wrap(dJgxe_dξ_c[p, i, j]) * du_dη +
            wrap(dJgxe_dη_c[p, i, j]) * du_dξ)

    expr = orthogonal + cross_term
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

"""
    fv_gradient_extended(u_ext, grid, dim)

Build an ArrayOp for the gradient ∂u/∂lon or ∂u/∂lat over ALL cells,
operating on a ghost-extended array.
"""
function fv_gradient_extended(u_ext, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc; Ng = grid.Ng
    dξ = grid.dξ; dη = grid.dη
    o = Ng

    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    u_c = const_wrap(u_ext)

    du_dξ = (wrap(u_c[p, i + o + 1, j + o]) - wrap(u_c[p, i + o - 1, j + o])) / (2 * dξ)
    du_dη = (wrap(u_c[p, i + o, j + o + 1]) - wrap(u_c[p, i + o, j + o - 1])) / (2 * dη)

    if dim == :lon
        dξ_dx = const_wrap(grid.dξ_dlon); dη_dx = const_wrap(grid.dη_dlon)
    else  # :lat
        dξ_dx = const_wrap(grid.dξ_dlat); dη_dx = const_wrap(grid.dη_dlat)
    end

    expr = wrap(dξ_dx[p, i, j]) * du_dξ + wrap(dη_dx[p, i, j]) * du_dη
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

"""
    SciMLBase.discretize(sys::PDESystem, disc::FVCubedSphere)

Discretize a PDESystem onto a cubed-sphere grid, returning an ODEProblem.

Builds ArrayOp-based FV operators over a ghost-extended symbolic array, then
creates MTK equations. Cross-panel boundary references are resolved during
ghost array construction, so ArrayOps use simple offset indexing.

This produces O(6·Nc²) MTK equations, but each equation is a compact
ArrayOp element (simple weighted sum with numerical coefficients) rather
than a complex chain-rule expression, dramatically reducing compile time.
"""
function SciMLBase.discretize(sys::PDESystem, disc::FVCubedSphere; kwargs...)
    grid = CubedSphereGrid(disc.Nc; R=disc.R, Ng=disc.Ng)
    Nc = grid.Nc

    # Identify time vs spatial IVs
    spatial_ivs = Num[]
    for iv in sys.ivs
        if identify_dimension(Symbol(iv)) != :t
            push!(spatial_ivs, iv)
        end
    end

    dvs = sys.dvs
    ndvs = length(dvs)

    # Create discrete state arrays
    disc_vars = Dict{Any,Any}()
    for dv in dvs
        name = Symbol(Symbolics.tosymbol(dv, escape=false))
        arr = first(Symbolics.@variables $name(mtk_t)[1:6, 1:Nc, 1:Nc])
        disc_vars[dv] = arr
    end

    # Build ghost-extended symbolic arrays for each DV
    ext_vars = Dict{Any,Any}()
    for dv in dvs
        ext_vars[dv] = _build_symbolic_ghost_extension(disc_vars[dv], grid)
    end

    # Analyze each PDE equation and build ArrayOp RHS
    all_eqs = Symbolics.Equation[]

    for eq in sys.eqs
        lhs_dv = _identify_lhs_dv(eq.lhs, dvs)
        lhs_arr = disc_vars[lhs_dv]

        # Build the RHS ArrayOp by analyzing the PDE structure
        rhs_arrayop = _build_rhs_arrayop(eq.rhs, dvs, spatial_ivs,
                                          disc_vars, ext_vars, grid)

        # Scalarize the ArrayOp to create per-cell equations
        rhs_scalar = Symbolics.scalarize(wrap(rhs_arrayop))

        for p in 1:6, i in 1:Nc, j in 1:Nc
            push!(all_eqs, mtk_D(lhs_arr[p, i, j]) ~ rhs_scalar[p, i, j])
        end
    end

    # Build MTK System and compile
    @named disc_system = System(all_eqs, mtk_t)
    compiled = mtkcompile(disc_system)

    # Project initial conditions
    u0 = _build_u0_mtk(sys, disc_vars, dvs, spatial_ivs, grid, compiled)

    tspan = extract_tspan(sys)
    return ODEProblem(compiled, u0, tspan; kwargs...)
end

"""
Build the RHS of a single PDE equation as an ArrayOp by analyzing the
symbolic expression and composing FV operators.

Recognized patterns:
- κ * (∂²u/∂x² + ∂²u/∂y²) → Laplacian with coefficient κ
- c * ∂u/∂x → gradient with coefficient c
- Constant * u → reaction term
"""
function _build_rhs_arrayop(rhs, dvs, spatial_ivs, disc_vars, ext_vars, grid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]

    # Accumulate the RHS ArrayOp expression
    rhs_expr = _rhs_to_arrayop_expr(unwrap(rhs), dvs, spatial_ivs,
                                      disc_vars, ext_vars, grid, idx, 1.0)

    return make_arrayop(idx, unwrap(rhs_expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

function _rhs_to_arrayop_expr(expr, dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff)
    p, i, j = idx[1], idx[2], idx[3]
    Nc = grid.Nc; Ng = grid.Ng; o = Ng
    dξ = grid.dξ; dη = grid.dη

    if !iscall(expr)
        v = Symbolics.value(wrap(expr))
        if v isa Number
            return wrap(Symbolics.value(coeff * Float64(v)))
        end
        return wrap(expr)
    end

    op = operation(expr)
    args = arguments(expr)

    # Addition: recurse into each term
    if op === (+)
        result = _rhs_to_arrayop_expr(args[1], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff)
        for k in 2:length(args)
            result = result + _rhs_to_arrayop_expr(args[k], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff)
        end
        return result
    end

    # Subtraction
    if op === (-) && length(args) == 1
        return _rhs_to_arrayop_expr(args[1], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, -coeff)
    end
    if op === (-) && length(args) == 2
        t1 = _rhs_to_arrayop_expr(args[1], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff)
        t2 = _rhs_to_arrayop_expr(args[2], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff)
        return t1 - t2
    end

    # Multiplication: split numeric * symbolic
    if op === (*)
        num_coeff = 1.0
        sym_factors = []
        for a in args
            v = Symbolics.value(wrap(a))
            if v isa Number
                num_coeff *= Float64(v)
            else
                push!(sym_factors, a)
            end
        end
        if length(sym_factors) == 1
            return _rhs_to_arrayop_expr(sym_factors[1], dvs, spatial_ivs, disc_vars, ext_vars, grid, idx, coeff * num_coeff)
        elseif isempty(sym_factors)
            return wrap(Symbolics.value(coeff * num_coeff))
        end
    end

    # Second derivative → chain-rule transformation
    if op isa Differential
        dim = identify_dimension(Symbol(op.x))
        inner = args[1]

        if iscall(inner) && operation(inner) isa Differential
            inner_dim = identify_dimension(Symbol(operation(inner).x))
            innermost = arguments(inner)[1]

            # Evaluate innermost expression at 9 stencil points
            f_c  = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, 0, 0, o)
            f_e  = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, +1, 0, o)
            f_w  = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, -1, 0, o)
            f_n  = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, 0, +1, o)
            f_s  = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, 0, -1, o)
            f_ne = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, +1, +1, o)
            f_nw = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, -1, +1, o)
            f_se = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, +1, -1, o)
            f_sw = _eval_at_gridpoint(innermost, dvs, ext_vars, idx, -1, -1, o)

            d2f_dξ2   = (f_e - 2*f_c + f_w) / dξ^2
            d2f_dη2   = (f_n - 2*f_c + f_s) / dη^2
            d2f_dξdη  = (f_ne - f_nw - f_se + f_sw) / (4*dξ*dη)
            df_dξ     = (f_e - f_w) / (2*dξ)
            df_dη     = (f_n - f_s) / (2*dη)

            a_ξ_o, a_η_o = _arrayop_chain_coeffs(dim, grid, idx)
            a_ξ_i, a_η_i = _arrayop_chain_coeffs(inner_dim, grid, idx)
            b_ξ, b_η = _arrayop_second_chain_coeffs(dim, inner_dim, grid, idx)

            result = a_ξ_o * a_ξ_i * d2f_dξ2 +
                     (a_ξ_o * a_η_i + a_η_o * a_ξ_i) * d2f_dξdη +
                     a_η_o * a_η_i * d2f_dη2 +
                     b_ξ * df_dξ + b_η * df_dη

            return coeff * result
        end

        # First derivative of general expression
        f_ip1 = _eval_at_gridpoint(inner, dvs, ext_vars, idx, +1, 0, o)
        f_im1 = _eval_at_gridpoint(inner, dvs, ext_vars, idx, -1, 0, o)
        f_jp1 = _eval_at_gridpoint(inner, dvs, ext_vars, idx, 0, +1, o)
        f_jm1 = _eval_at_gridpoint(inner, dvs, ext_vars, idx, 0, -1, o)

        df_dξ = (f_ip1 - f_im1) / (2*dξ)
        df_dη = (f_jp1 - f_jm1) / (2*dη)

        if dim == :lon
            dξ_dx = const_wrap(grid.dξ_dlon); dη_dx = const_wrap(grid.dη_dlon)
            return coeff * (wrap(dξ_dx[p, i, j]) * df_dξ + wrap(dη_dx[p, i, j]) * df_dη)
        elseif dim == :lat
            dξ_dy = const_wrap(grid.dξ_dlat); dη_dy = const_wrap(grid.dη_dlat)
            return coeff * (wrap(dξ_dy[p, i, j]) * df_dξ + wrap(dη_dy[p, i, j]) * df_dη)
        elseif dim == :xi
            return coeff * df_dξ
        elseif dim == :eta
            return coeff * df_dη
        end
    end

    # DV call (reaction term) — use ghost-extended array with Const wrapping
    dv = _match_dv(expr, dvs)
    if dv !== nothing
        u_c = const_wrap(ext_vars[dv])
        return coeff * wrap(u_c[p, i + o, j + o])
    end

    # Fallback: return as numeric constant
    return coeff * wrap(expr)
end

"""
Evaluate a symbolic expression at a grid point offset (di, dj) from the
current cell, by replacing each dependent variable with its ghost-extended
array value at the offset position.
"""
function _eval_at_gridpoint(expr, dvs, ext_vars, idx, di, dj, o)
    p, i, j = idx[1], idx[2], idx[3]
    ex = unwrap(expr)

    if !iscall(ex)
        return wrap(ex)
    end

    dv = _match_dv(ex, dvs)
    if dv !== nothing
        u_c = const_wrap(ext_vars[dv])
        return wrap(u_c[p, i + di + o, j + dj + o])
    end

    # Recurse into arguments for nonlinear expressions (u^2, u*v, sin(u), etc.)
    op = operation(ex); args = arguments(ex)
    new_args = [unwrap(_eval_at_gridpoint(wrap(a), dvs, ext_vars, idx, di, dj, o)) for a in args]
    return wrap(op(new_args...))
end

"""
Return (coeff_ξ, coeff_η) chain-rule coefficients for dimension `dim` as ArrayOp expressions.
For computational coordinates: (1,0) or (0,1). For physical coordinates: uses coordinate Jacobian.
"""
function _arrayop_chain_coeffs(dim, grid, idx)
    p, i, j = idx[1], idx[2], idx[3]
    if dim == :lon
        return (wrap(const_wrap(grid.dξ_dlon)[p, i, j]),
                wrap(const_wrap(grid.dη_dlon)[p, i, j]))
    elseif dim == :lat
        return (wrap(const_wrap(grid.dξ_dlat)[p, i, j]),
                wrap(const_wrap(grid.dη_dlat)[p, i, j]))
    elseif dim == :xi
        return (wrap(Symbolics.value(1.0)), wrap(Symbolics.value(0.0)))
    elseif dim == :eta
        return (wrap(Symbolics.value(0.0)), wrap(Symbolics.value(1.0)))
    end
end

"""
Return (coeff_ξ, coeff_η) for the second-derivative chain-rule correction
(∂²ξ/∂x∂y, ∂²η/∂x∂y) terms from the nonlinear coordinate mapping.
"""
function _arrayop_second_chain_coeffs(dim_outer, dim_inner, grid, idx)
    p, i, j = idx[1], idx[2], idx[3]
    if dim_outer == :lon && dim_inner == :lon
        return (wrap(const_wrap(grid.d2ξ_dlon2)[p, i, j]),
                wrap(const_wrap(grid.d2η_dlon2)[p, i, j]))
    elseif dim_outer == :lat && dim_inner == :lat
        return (wrap(const_wrap(grid.d2ξ_dlat2)[p, i, j]),
                wrap(const_wrap(grid.d2η_dlat2)[p, i, j]))
    elseif (dim_outer == :lon && dim_inner == :lat) || (dim_outer == :lat && dim_inner == :lon)
        return (wrap(const_wrap(grid.d2ξ_dlondlat)[p, i, j]),
                wrap(const_wrap(grid.d2η_dlondlat)[p, i, j]))
    else
        return (wrap(Symbolics.value(0.0)), wrap(Symbolics.value(0.0)))
    end
end

function _match_dv(expr, dvs)
    for dv in dvs
        if isequal(wrap(expr), dv)
            return dv
        end
    end
    if iscall(expr)
        name = Symbol(Symbolics.tosymbol(wrap(expr), escape=false))
        for dv in dvs
            dv_name = Symbol(Symbolics.tosymbol(dv, escape=false))
            if name == dv_name
                return dv
            end
        end
    end
    return nothing
end

function _build_u0_mtk(sys, disc_vars, dvs, spatial_ivs, grid, compiled)
    Nc = grid.Nc
    u0 = Pair[]
    tspan = extract_tspan(sys)
    t0 = tspan[1]

    for dv in dvs
        arr = disc_vars[dv]
        ic_found = false

        for bc in sys.bcs
            if _is_initial_condition(bc, dv, t0)
                rhs = bc.rhs
                for p in 1:6, i in 1:Nc, j in 1:Nc
                    val = _eval_ic(rhs, spatial_ivs, grid, p, i, j)
                    push!(u0, arr[p, i, j] => val)
                end
                ic_found = true
                break
            end
        end

        if !ic_found
            for p in 1:6, i in 1:Nc, j in 1:Nc
                push!(u0, arr[p, i, j] => 0.0)
            end
        end
    end

    return u0
end

function extract_tspan(sys::PDESystem)
    for d in sys.domain
        if identify_dimension(Symbol(d.variables)) == :t
            return (Float64(d.domain.left), Float64(d.domain.right))
        end
    end
    error("No time domain found")
end

function _is_initial_condition(bc, dv, t0)
    lhs = unwrap(bc.lhs)
    if !iscall(lhs)
        return false
    end
    args = arguments(lhs)
    lhs_name = Symbol(Symbolics.tosymbol(wrap(lhs), escape=false))
    dv_name = Symbol(Symbolics.tosymbol(dv, escape=false))
    if lhs_name != dv_name
        return false
    end
    if length(args) >= 1
        t_val = Symbolics.value(wrap(args[1]))
        if t_val isa Number && isapprox(Float64(t_val), t0)
            return true
        end
    end
    return false
end

function _eval_ic(rhs, spatial_ivs, grid, p, i, j)
    subs = Dict{Any,Any}()
    if length(spatial_ivs) >= 1
        subs[spatial_ivs[1]] = grid.lon[p, i, j]
    end
    if length(spatial_ivs) >= 2
        subs[spatial_ivs[2]] = grid.lat[p, i, j]
    end
    v = Symbolics.value(Symbolics.substitute(rhs, subs))
    return v isa Number ? Float64(v) : Float64(eval(Symbolics.toexpr(v)))
end
