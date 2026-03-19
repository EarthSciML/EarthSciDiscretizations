"""
1D flux-form transport operators.

Provides both the original Lax-Friedrichs ArrayOp operator (for symbolic use)
and a PPM-based numerical operator (for accurate transport).
"""

"""
    flux_1d(q, courant, grid, dim)

ArrayOp for 1D Lax-Friedrichs transport tendency. Index (p,i,j) maps to
physical cell (p, i+1, j) for `:xi` or (p, i, j+1) for `:eta`.

Note: This is a first-order scheme retained for backward compatibility and
symbolic testing. For accurate transport, use `flux_1d_ppm!`.
"""
function flux_1d(q, courant, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc; idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    q_c = const_wrap(unwrap(q)); c_c = const_wrap(unwrap(courant))
    if dim == :xi
        dξ = grid.dξ
        F_e_a = wrap(c_c[p, i + 2, j]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 2, j])) / 2
        F_e_d = abs(wrap(c_c[p, i + 2, j])) * (wrap(q_c[p, i + 2, j]) - wrap(q_c[p, i + 1, j])) / 2
        F_w_a = wrap(c_c[p, i + 1, j]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i + 1, j])) / 2
        F_w_d = abs(wrap(c_c[p, i + 1, j])) * (wrap(q_c[p, i + 1, j]) - wrap(q_c[p, i, j])) / 2
        expr = -((F_e_a - F_e_d) - (F_w_a - F_w_d)) / dξ
        ranges = Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:Nc)
    else
        dη = grid.dη
        F_n_a = wrap(c_c[p, i, j + 2]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i, j + 2])) / 2
        F_n_d = abs(wrap(c_c[p, i, j + 2])) * (wrap(q_c[p, i, j + 2]) - wrap(q_c[p, i, j + 1])) / 2
        F_s_a = wrap(c_c[p, i, j + 1]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i, j + 1])) / 2
        F_s_d = abs(wrap(c_c[p, i, j + 1])) * (wrap(q_c[p, i, j + 1]) - wrap(q_c[p, i, j])) / 2
        expr = -((F_n_a - F_n_d) - (F_s_a - F_s_d)) / dη
        ranges = Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc - 2))
    end
    return make_arrayop(idx, unwrap(expr), ranges)
end

"""
    flux_1d_ppm!(tendency, q, vel, grid, dim, dt)

PPM-based 1D flux-form transport. Computes the tendency for `q` transported
by velocity `vel` in dimension `dim` (:xi or :eta).

Uses PPM reconstruction with the Colella-Woodward (1984) limiter for
high-order, monotone transport. The flux through each interface is computed
by integrating the parabolic profile over the swept volume.

Arguments:
- `tendency`: output array (6, Nc, Nc), modified in-place
- `q`: scalar field (6, Nc, Nc)
- `vel`: velocity field at cell edges (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `grid`: CubedSphereGrid
- `dim`: :xi or :eta
- `dt`: time step
"""
function flux_1d_ppm!(tendency, q, vel, grid::CubedSphereGrid, dim::Symbol, dt)
    Nc = grid.Nc

    # Extend q with ghost cells for stencil access at boundaries
    q_ext = extend_with_ghosts(q, grid)
    Ng = grid.Ng

    if dim == :xi
        # Compute interface fluxes in the ξ-direction
        # Interface at edge i sits between cell i-1 and cell i
        flux = zeros(6, Nc + 1, Nc)
        for p in 1:6, j in 1:Nc
            for i in 1:Nc+1
                # Cells adjacent to this interface (in extended array)
                ie = i + Ng  # extended array index for cell i (or edge i)
                je = j + Ng

                # PPM reconstruction for the upwind cell
                # We need cell values: ie-2, ie-1, ie, ie+1, ie+2 in extended array
                qim2 = q_ext[p, ie - 2, je]; qim1 = q_ext[p, ie - 1, je]
                qi   = q_ext[p, ie, je];     qip1 = q_ext[p, ie + 1, je]
                qip2 = q_ext[p, ie + 2, je]

                # Interface value at i-1/2 (between cell i-1 and cell i)
                qi_half = (7.0 / 12.0) * (qim1 + qi) - (1.0 / 12.0) * (qim2 + qip1)

                # Compute Courant number at this interface
                c = vel[p, i, j] * dt / (grid.R * grid.dξ)

                if c >= 0
                    # Upwind cell is i-1 (to the left)
                    # Need left/right edges of cell i-1
                    ql_left = (7.0 / 12.0) * (qim2 + qim1) - (1.0 / 12.0) * (q_ext[p, ie - 3, je] + qi)
                    qr_left = qi_half
                    ql_left, qr_left = _ppm_limit_cw84(ql_left, qr_left, qim1)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_left, qr_left, qim1, c)
                else
                    # Upwind cell is i (to the right)
                    ql_right = qi_half
                    qr_right = (7.0 / 12.0) * (qi + qip1) - (1.0 / 12.0) * (qim1 + qip2)
                    ql_right, qr_right = _ppm_limit_cw84(ql_right, qr_right, qi)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_right, qr_right, qi, c)
                end
            end
        end

        # Compute tendency: -(F_{i+1} - F_i) / (A * dξ)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(flux[p, i + 1, j] * grid.dx[p, i + 1, j] -
                                   flux[p, i, j] * grid.dx[p, i, j]) / grid.area[p, i, j]
        end
    else  # dim == :eta
        flux = zeros(6, Nc, Nc + 1)
        for p in 1:6, i in 1:Nc
            for j in 1:Nc+1
                ie = i + Ng; je = j + Ng

                qjm2 = q_ext[p, ie, je - 2]; qjm1 = q_ext[p, ie, je - 1]
                qj   = q_ext[p, ie, je];     qjp1 = q_ext[p, ie, je + 1]
                qjp2 = q_ext[p, ie, je + 2]

                qj_half = (7.0 / 12.0) * (qjm1 + qj) - (1.0 / 12.0) * (qjm2 + qjp1)

                c = vel[p, i, j] * dt / (grid.R * grid.dη)

                if c >= 0
                    ql_left = (7.0 / 12.0) * (qjm2 + qjm1) - (1.0 / 12.0) * (q_ext[p, ie, je - 3] + qj)
                    qr_left = qj_half
                    ql_left, qr_left = _ppm_limit_cw84(ql_left, qr_left, qjm1)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_left, qr_left, qjm1, c)
                else
                    ql_right = qj_half
                    qr_right = (7.0 / 12.0) * (qj + qjp1) - (1.0 / 12.0) * (qjm1 + qjp2)
                    ql_right, qr_right = _ppm_limit_cw84(ql_right, qr_right, qj)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_right, qr_right, qj, c)
                end
            end
        end

        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(flux[p, i, j + 1] * grid.dy[p, i, j + 1] -
                                   flux[p, i, j] * grid.dy[p, i, j]) / grid.area[p, i, j]
        end
    end

    return tendency
end
