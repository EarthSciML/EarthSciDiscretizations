"""
Two-sided PPM interpolation at cube edges (FV3 Section 6.5).

At cube edges, the standard 4th-order PPM interface interpolation formula
breaks down because the coordinate system is discontinuous across panel
boundaries. FV3 uses a two-sided extrapolation that accounts for varying
grid-cell widths:

    a_E = (1/2) · [(2dx₁ + dx₂)q₁ - dx₁ q₂] / (dx₁ + dx₂)
        + (1/2) · [(2dx₀ + dx₋₁)q₀ - dx₀ q₋₁] / (dx₀ + dx₋₁)

where indices ≤ 0 are on the neighboring panel, accessed through ghost cells.
This produces a smooth transition across the discontinuous coordinate system
at panel boundaries.

For monotonic schemes, an additional constraint bounds the edge value within
the range of the 4 cell-mean values (FV3 Eq. 6.6).

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Eq. 6.5-6.6.
"""

"""
    ppm_edge_value_twosided(q0, q1, q2, q_m1, dx0, dx1, dx2, dx_m1)

Compute the PPM interface value at a cube edge using two-sided extrapolation.

Arguments:
- `q0`, `q_m1`: cell-mean values on the "minus" side (neighbor panel, through ghost cells)
- `q1`, `q2`: cell-mean values on the "plus" side (local panel)
- `dx0`, `dx_m1`: grid cell widths on the minus side
- `dx1`, `dx2`: grid cell widths on the plus side

The interface is between q0 (minus side, index 0) and q1 (plus side, index 1).

Returns the two-sided interface value a_E.
"""
function ppm_edge_value_twosided(q0, q1, q2, q_m1, dx0, dx1, dx2, dx_m1)
    # Right-side extrapolation: from cells 1 and 2 toward the edge
    right = (2 * dx1 + dx2) * q1 - dx1 * q2
    right /= (dx1 + dx2)

    # Left-side extrapolation: from cells 0 and -1 toward the edge
    left = (2 * dx0 + dx_m1) * q0 - dx0 * q_m1
    left /= (dx0 + dx_m1)

    # Average the two one-sided estimates
    a_E = 0.5 * (right + left)

    return a_E
end

"""
    ppm_edge_value_twosided_limited(q0, q1, q2, q_m1, dx0, dx1, dx2, dx_m1)

Two-sided PPM edge value with monotonicity constraint (FV3 Eq. 6.6).

After computing the two-sided estimate, clamp the result to lie within
the range [min(q0, q1, q2, q_m1), max(q0, q1, q2, q_m1)] to prevent
new extrema from forming at cube edges.
"""
function ppm_edge_value_twosided_limited(q0, q1, q2, q_m1, dx0, dx1, dx2, dx_m1)
    a_E = ppm_edge_value_twosided(q0, q1, q2, q_m1, dx0, dx1, dx2, dx_m1)

    # Monotonicity constraint: bound within the range of the 4 contributing cells
    qmin = min(q0, q1, q2, q_m1)
    qmax = max(q0, q1, q2, q_m1)
    a_E = clamp(a_E, qmin, qmax)

    return a_E
end

"""
    flux_1d_ppm_twosided!(tendency, q, vel, grid, dim, dt)

PPM-based 1D flux-form transport with two-sided interpolation at cube edges.

This is an enhanced version of `flux_1d_ppm!` that uses the FV3 two-sided
extrapolation (Eq. 6.5) at panel boundary interfaces, where the standard
4th-order PPM formula is inaccurate due to coordinate discontinuity.

At interior interfaces (≥ 3 cells from panel edges), the standard 4th-order
formula is used. At interfaces within 2 cells of a panel boundary, the
two-sided formula with explicit grid-cell widths is used instead.
"""
function flux_1d_ppm_twosided!(tendency, q, vel, grid::CubedSphereGrid, dim::Symbol, dt)
    Nc = grid.Nc
    Ng = grid.Ng

    # Extend q with ghost cells
    q_ext = extend_with_ghosts(q, grid)

    if dim == :xi
        flux = zeros(6, Nc + 1, Nc)

        # Precompute cell widths including ghost cells
        # Cell width in ξ at cell center: approximate using edge-to-edge distance
        dxc = zeros(6, Nc + 2 * Ng)
        for p in 1:6
            for i in 1:Nc
                dxc[p, i + Ng] = grid.dξ  # uniform in computational space
            end
            # Ghost cell widths: same as interior (gnomonic grid has uniform dξ)
            for g in 1:Ng
                dxc[p, g] = grid.dξ
                dxc[p, Nc + Ng + g] = grid.dξ
            end
        end

        for p in 1:6, j in 1:Nc
            for i in 1:(Nc + 1)
                ie = i + Ng  # extended index for interface between cell ie-1 and ie

                # Determine if this interface is near a panel boundary
                near_boundary = (i <= 2) || (i >= Nc)

                # Get stencil values from extended array
                qim2 = q_ext[p, ie - 2, j + Ng]
                qim1 = q_ext[p, ie - 1, j + Ng]
                qi = q_ext[p, ie, j + Ng]
                qip1 = q_ext[p, ie + 1, j + Ng]

                if near_boundary
                    # Two-sided extrapolation at/near panel boundaries
                    # Interface between cell ie-1 (q0=qim1) and cell ie (q1=qi)
                    qi_half = ppm_edge_value_twosided_limited(
                        qim1, qi, qip1, qim2,
                        dxc[p, ie - 1], dxc[p, ie], dxc[p, ie + 1], dxc[p, ie - 2]
                    )
                else
                    # Standard 4th-order formula for interior
                    qi_half = (7.0 / 12.0) * (qim1 + qi) - (1.0 / 12.0) * (qim2 + qip1)
                end

                # Compute Courant number
                c = _get_courant_xi(vel, dt, grid, p, i, j)

                if c >= 0
                    # Upwind cell is ie-1
                    qim3 = q_ext[p, ie - 3, j + Ng]
                    if near_boundary || i <= 3
                        ql_left = ppm_edge_value_twosided_limited(
                            qim2, qim1, qi, qim3,
                            dxc[p, ie - 2], dxc[p, ie - 1], dxc[p, ie], dxc[p, ie - 3]
                        )
                    else
                        ql_left = (7.0 / 12.0) * (qim2 + qim1) - (1.0 / 12.0) * (qim3 + qi)
                    end
                    qr_left = qi_half
                    ql_left, qr_left = _ppm_limit_cw84(ql_left, qr_left, qim1)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_left, qr_left, qim1, c)
                else
                    # Upwind cell is ie
                    qip2 = q_ext[p, ie + 2, j + Ng]
                    ql_right = qi_half
                    if near_boundary || i >= Nc - 1
                        qr_right = ppm_edge_value_twosided_limited(
                            qi, qip1, qip2, qim1,
                            dxc[p, ie], dxc[p, ie + 1], dxc[p, ie + 2], dxc[p, ie - 1]
                        )
                    else
                        qr_right = (7.0 / 12.0) * (qi + qip1) - (1.0 / 12.0) * (qim1 + qip2)
                    end
                    ql_right, qr_right = _ppm_limit_cw84(ql_right, qr_right, qi)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_right, qr_right, qi, c)
                end
            end
        end

        _match_boundary_fluxes_xi!(flux, grid)

        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(
                flux[p, i + 1, j] * grid.dx[p, i + 1, j] -
                    flux[p, i, j] * grid.dx[p, i, j]
            ) / grid.area[p, i, j]
        end

    else  # dim == :eta
        flux = zeros(6, Nc, Nc + 1)

        dxc = zeros(6, Nc + 2 * Ng)
        for p in 1:6
            for j in 1:Nc
                dxc[p, j + Ng] = grid.dη
            end
            for g in 1:Ng
                dxc[p, g] = grid.dη
                dxc[p, Nc + Ng + g] = grid.dη
            end
        end

        for p in 1:6, i in 1:Nc
            for j in 1:(Nc + 1)
                je = j + Ng

                near_boundary = (j <= 2) || (j >= Nc)

                qjm2 = q_ext[p, i + Ng, je - 2]
                qjm1 = q_ext[p, i + Ng, je - 1]
                qj = q_ext[p, i + Ng, je]
                qjp1 = q_ext[p, i + Ng, je + 1]

                if near_boundary
                    qj_half = ppm_edge_value_twosided_limited(
                        qjm1, qj, qjp1, qjm2,
                        dxc[p, je - 1], dxc[p, je], dxc[p, je + 1], dxc[p, je - 2]
                    )
                else
                    qj_half = (7.0 / 12.0) * (qjm1 + qj) - (1.0 / 12.0) * (qjm2 + qjp1)
                end

                c = _get_courant_eta(vel, dt, grid, p, i, j)

                if c >= 0
                    qjm3 = q_ext[p, i + Ng, je - 3]
                    if near_boundary || j <= 3
                        ql_left = ppm_edge_value_twosided_limited(
                            qjm2, qjm1, qj, qjm3,
                            dxc[p, je - 2], dxc[p, je - 1], dxc[p, je], dxc[p, je - 3]
                        )
                    else
                        ql_left = (7.0 / 12.0) * (qjm2 + qjm1) - (1.0 / 12.0) * (qjm3 + qj)
                    end
                    qr_left = qj_half
                    ql_left, qr_left = _ppm_limit_cw84(ql_left, qr_left, qjm1)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_left, qr_left, qjm1, c)
                else
                    qjp2 = q_ext[p, i + Ng, je + 2]
                    ql_right = qj_half
                    if near_boundary || j >= Nc - 1
                        qr_right = ppm_edge_value_twosided_limited(
                            qj, qjp1, qjp2, qjm1,
                            dxc[p, je], dxc[p, je + 1], dxc[p, je + 2], dxc[p, je - 1]
                        )
                    else
                        qr_right = (7.0 / 12.0) * (qj + qjp1) - (1.0 / 12.0) * (qjm1 + qjp2)
                    end
                    ql_right, qr_right = _ppm_limit_cw84(ql_right, qr_right, qj)
                    flux[p, i, j] = vel[p, i, j] * ppm_flux_integral(ql_right, qr_right, qj, c)
                end
            end
        end

        _match_boundary_fluxes_eta!(flux, grid)

        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(
                flux[p, i, j + 1] * grid.dy[p, i, j + 1] -
                    flux[p, i, j] * grid.dy[p, i, j]
            ) / grid.area[p, i, j]
        end
    end

    return tendency
end
