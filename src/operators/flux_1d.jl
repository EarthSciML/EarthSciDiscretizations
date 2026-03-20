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
    _get_courant_xi(vel, dt, grid, p, i, j)

Compute the Courant number for the ξ-direction interface at (p, i, j),
using proper boundary distances at panel edges.
"""
function _get_courant_xi(vel, dt, grid::CubedSphereGrid, p, i, j)
    Nc = grid.Nc
    if i == 1
        dist = grid.dist_xi_bnd[p, 1, j]
    elseif i == Nc + 1
        dist = grid.dist_xi_bnd[p, 2, j]
    else
        dist = grid.dist_xi[p, i - 1, j]
    end
    return vel[p, i, j] * dt / dist
end

"""
    _get_courant_eta(vel, dt, grid, p, i, j)

Compute the Courant number for the η-direction interface at (p, i, j),
using proper boundary distances at panel edges.
"""
function _get_courant_eta(vel, dt, grid::CubedSphereGrid, p, i, j)
    Nc = grid.Nc
    if j == 1
        dist = grid.dist_eta_bnd[p, i, 1]
    elseif j == Nc + 1
        dist = grid.dist_eta_bnd[p, i, 2]
    else
        dist = grid.dist_eta[p, i, j - 1]
    end
    return vel[p, i, j] * dt / dist
end

"""
    _compute_ppm_fluxes!(flux, q, vel, grid, dim, dt)

Compute PPM fluxes at cell interfaces without converting to tendencies.
Returns the flux array directly for external flux matching.

Uses PPM reconstruction with the Colella-Woodward (1984) limiter for
high-order, monotone transport. The flux through each interface is computed
by integrating the parabolic profile over the swept volume.

At panel boundaries, same-direction fluxes (ξ-to-ξ or η-to-η) are matched
between adjacent panels to ensure conservation. For rotated panel boundary
connections, use `_match_rotated_boundary_fluxes!` after computing both
ξ and η fluxes.

Arguments:
- `flux`: output flux array, (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `q`: scalar field (6, Nc, Nc)
- `vel`: velocity field at cell edges (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `grid`: CubedSphereGrid
- `dim`: :xi or :eta
- `dt`: time step
"""
function _compute_ppm_fluxes!(flux, q, vel, grid::CubedSphereGrid, dim::Symbol, dt)
    Nc = grid.Nc

    # Extend q with ghost cells for stencil access at boundaries
    q_ext = extend_with_ghosts(q, grid)
    Ng = grid.Ng

    if dim == :xi
        # Compute interface fluxes in the ξ-direction
        # Interface at edge i sits between cell i-1 and cell i
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

                # Compute Courant number using proper boundary distances
                c = _get_courant_xi(vel, dt, grid, p, i, j)

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

        # Match boundary fluxes for same-direction panel connections
        _match_boundary_fluxes_xi!(flux, grid)
    else  # dim == :eta
        for p in 1:6, i in 1:Nc
            for j in 1:Nc+1
                ie = i + Ng; je = j + Ng

                qjm2 = q_ext[p, ie, je - 2]; qjm1 = q_ext[p, ie, je - 1]
                qj   = q_ext[p, ie, je];     qjp1 = q_ext[p, ie, je + 1]
                qjp2 = q_ext[p, ie, je + 2]

                qj_half = (7.0 / 12.0) * (qjm1 + qj) - (1.0 / 12.0) * (qjm2 + qjp1)

                # Compute Courant number using proper boundary distances
                c = _get_courant_eta(vel, dt, grid, p, i, j)

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

        # Match boundary fluxes for same-direction panel connections
        _match_boundary_fluxes_eta!(flux, grid)
    end

    return flux
end

"""
    _flux_to_tendency!(tendency, flux, grid, dim)

Convert a flux array to tendency by computing the divergence:
    tendency[p,i,j] = -(F_{i+1} * edge_length_{i+1} - F_i * edge_length_i) / area

Arguments:
- `tendency`: output array (6, Nc, Nc), modified in-place
- `flux`: flux array, (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `grid`: CubedSphereGrid
- `dim`: :xi or :eta
"""
function _flux_to_tendency!(tendency, flux, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(flux[p, i + 1, j] * grid.dx[p, i + 1, j] -
                                   flux[p, i, j] * grid.dx[p, i, j]) / grid.area[p, i, j]
        end
    else
        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(flux[p, i, j + 1] * grid.dy[p, i, j + 1] -
                                   flux[p, i, j] * grid.dy[p, i, j]) / grid.area[p, i, j]
        end
    end
    return tendency
end

"""
    flux_1d_ppm!(tendency, q, vel, grid, dim, dt)

PPM-based 1D flux-form transport. Computes the tendency for `q` transported
by velocity `vel` in dimension `dim` (:xi or :eta).

Uses PPM reconstruction with the Colella-Woodward (1984) limiter for
high-order, monotone transport. The flux through each interface is computed
by integrating the parabolic profile over the swept volume.

At panel boundaries, fluxes are matched between adjacent panels to ensure
exact conservation: for same-direction connections, the flux is averaged
between both panels' independently computed values.

Arguments:
- `tendency`: output array (6, Nc, Nc), modified in-place
- `q`: scalar field (6, Nc, Nc)
- `vel`: velocity field at cell edges (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `grid`: CubedSphereGrid
- `dim`: :xi or :eta
- `dt`: time step

!!! note "Conservation at panel boundaries"
    Same-direction boundary fluxes (ξ-to-ξ, η-to-η) are matched within this
    function. For rotated boundary connections (ξ-to-η or η-to-ξ), use
    `_compute_ppm_fluxes!` to get the raw flux arrays, then call
    `_match_rotated_boundary_fluxes!` before converting to tendencies.
"""
function flux_1d_ppm!(tendency, q, vel, grid::CubedSphereGrid, dim::Symbol, dt)
    Nc = grid.Nc
    if dim == :xi
        flux = zeros(6, Nc + 1, Nc)
    else
        flux = zeros(6, Nc, Nc + 1)
    end
    _compute_ppm_fluxes!(flux, q, vel, grid, dim, dt)
    _flux_to_tendency!(tendency, flux, grid, dim)
    return tendency
end

"""
    _match_boundary_fluxes_xi!(flux, grid)

Match ξ-direction fluxes at panel boundaries for same-direction connections.
For each shared edge where both panels have ξ-fluxes, average the two
independently computed fluxes to ensure exact conservation.
"""
function _match_boundary_fluxes_xi!(flux, grid::CubedSphereGrid)
    Nc = grid.Nc
    for p in 1:6
        # East boundary: flux[p, Nc+1, j] is the east-edge flux of panel p
        nb = PANEL_CONNECTIVITY[p][East]
        nb_edge = nb.neighbor_edge
        # Only match if neighbor's edge is West (same-direction ξ-ξ connection)
        if nb_edge == West
            for j in 1:Nc
                j_nb = nb.reverse_index ? (Nc + 1 - j) : j
                # Panel p's east flux = flux[p, Nc+1, j]
                # Neighbor's west flux = flux[nb.neighbor_panel, 1, j_nb]
                # The fluxes should be equal in magnitude (outgoing from p = incoming to neighbor)
                # Average to enforce exact conservation
                avg = 0.5 * (flux[p, Nc + 1, j] + flux[nb.neighbor_panel, 1, j_nb])
                flux[p, Nc + 1, j] = avg
                flux[nb.neighbor_panel, 1, j_nb] = avg
            end
        end
    end
end

"""
    _match_boundary_fluxes_eta!(flux, grid)

Match η-direction fluxes at panel boundaries for same-direction connections.
"""
function _match_boundary_fluxes_eta!(flux, grid::CubedSphereGrid)
    Nc = grid.Nc
    for p in 1:6
        # North boundary: flux[p, i, Nc+1]
        nb = PANEL_CONNECTIVITY[p][North]
        nb_edge = nb.neighbor_edge
        # Only match if neighbor's edge is South (same-direction η-η connection)
        if nb_edge == South
            for i in 1:Nc
                i_nb = nb.reverse_index ? (Nc + 1 - i) : i
                avg = 0.5 * (flux[p, i, Nc + 1] + flux[nb.neighbor_panel, i_nb, 1])
                flux[p, i, Nc + 1] = avg
                flux[nb.neighbor_panel, i_nb, 1] = avg
            end
        end
    end
end

"""
    _get_boundary_mass_flux(flux_xi, flux_eta, grid, p, edge, k)

Get the outward mass flux through the k-th cell interface on `edge` of panel `p`.
Positive return value means flow OUT of panel p.

- For East: mass_flux = flux_xi[p, Nc+1, k] * dx[p, Nc+1, k]
- For West: mass_flux = -flux_xi[p, 1, k] * dx[p, 1, k]  (positive xi = inward)
- For North: mass_flux = flux_eta[p, k, Nc+1] * dy[p, k, Nc+1]
- For South: mass_flux = -flux_eta[p, k, 1] * dy[p, k, 1]  (positive eta = inward)
"""
function _get_boundary_mass_flux(flux_xi, flux_eta, grid::CubedSphereGrid, p, edge::EdgeDirection, k)
    Nc = grid.Nc
    if edge == East
        return flux_xi[p, Nc + 1, k] * grid.dx[p, Nc + 1, k]
    elseif edge == West
        return -flux_xi[p, 1, k] * grid.dx[p, 1, k]
    elseif edge == North
        return flux_eta[p, k, Nc + 1] * grid.dy[p, k, Nc + 1]
    else  # South
        return -flux_eta[p, k, 1] * grid.dy[p, k, 1]
    end
end

"""
    _set_boundary_mass_flux!(flux_xi, flux_eta, grid, p, edge, k, mass_flux_out)

Set the boundary flux value such that the outward mass flux equals `mass_flux_out`.

Inverts the sign convention from `_get_boundary_mass_flux`:
- East: flux_xi[p, Nc+1, k] = mass_flux_out / dx
- West: flux_xi[p, 1, k] = -mass_flux_out / dx
- North: flux_eta[p, k, Nc+1] = mass_flux_out / dy
- South: flux_eta[p, k, 1] = -mass_flux_out / dy
"""
function _set_boundary_mass_flux!(flux_xi, flux_eta, grid::CubedSphereGrid, p, edge::EdgeDirection, k, mass_flux_out)
    Nc = grid.Nc
    if edge == East
        flux_xi[p, Nc + 1, k] = mass_flux_out / grid.dx[p, Nc + 1, k]
    elseif edge == West
        flux_xi[p, 1, k] = -mass_flux_out / grid.dx[p, 1, k]
    elseif edge == North
        flux_eta[p, k, Nc + 1] = mass_flux_out / grid.dy[p, k, Nc + 1]
    else  # South
        flux_eta[p, k, 1] = -mass_flux_out / grid.dy[p, k, 1]
    end
end

"""
    _match_rotated_boundary_fluxes!(flux_xi, flux_eta, grid)

Match fluxes at panel boundaries that are NOT handled by the same-direction
matchers `_match_boundary_fluxes_xi!` and `_match_boundary_fluxes_eta!`.

This includes:
- Rotated connections where a ξ-edge of one panel connects to an η-edge
  of its neighbor (or vice versa)
- Anti-parallel same-direction connections (e.g., North-to-North, South-to-South)

For conservation, the mass flux leaving one panel through a shared edge must
equal the mass flux entering the adjacent panel through that same edge.
Each shared edge is processed exactly once by using panel ordering to avoid
double-processing.

The outward mass flux convention is:
- East:  +flux_xi[p, Nc+1, k] * dx   (positive ξ-flux = outward)
- West:  -flux_xi[p, 1, k] * dx      (positive ξ-flux = inward, so negate)
- North: +flux_eta[p, k, Nc+1] * dy  (positive η-flux = outward)
- South: -flux_eta[p, k, 1] * dy     (positive η-flux = inward, so negate)

At a shared edge: M_out_p + M_out_q = 0 (outflow from p = inflow to q).
We average: avg = 0.5 * (M_out_p - M_out_q), then set M_out_p = avg, M_out_q = -avg.
"""
function _match_rotated_boundary_fluxes!(flux_xi, flux_eta, grid::CubedSphereGrid)
    Nc = grid.Nc

    # Track which physical edges have been processed to avoid double-processing.
    # Key: (min_panel, max_panel, min_edge, max_edge) or similar unique identifier.
    processed = Set{Tuple{Int,EdgeDirection,Int,EdgeDirection}}()

    for p in 1:6
        for edge_p in (East, West, North, South)
            nb = PANEL_CONNECTIVITY[p][edge_p]
            nb_panel = nb.neighbor_panel
            nb_edge = nb.neighbor_edge

            # Skip same-direction connections already handled:
            # East→West handled by _match_boundary_fluxes_xi!
            # West→East is the reverse side of an East→West (also handled)
            # North→South handled by _match_boundary_fluxes_eta!
            # South→North is the reverse side of a North→South (also handled)
            if (edge_p == East && nb_edge == West) ||
               (edge_p == West && nb_edge == East) ||
               (edge_p == North && nb_edge == South) ||
               (edge_p == South && nb_edge == North)
                continue
            end

            # Create a canonical key for this physical edge to avoid double-processing
            edge_key = if p < nb_panel
                (p, edge_p, nb_panel, nb_edge)
            elseif p > nb_panel
                (nb_panel, nb_edge, p, edge_p)
            else
                # Same panel (shouldn't happen in standard cube, but handle it)
                if edge_p <= nb_edge
                    (p, edge_p, nb_panel, nb_edge)
                else
                    (nb_panel, nb_edge, p, edge_p)
                end
            end

            edge_key in processed && continue
            push!(processed, edge_key)

            # Process each cell along this shared edge
            for k in 1:Nc
                k_nb = nb.reverse_index ? (Nc + 1 - k) : k

                # Get outward mass fluxes from each panel's perspective
                M_out_p = _get_boundary_mass_flux(flux_xi, flux_eta, grid, p, edge_p, k)
                M_out_q = _get_boundary_mass_flux(flux_xi, flux_eta, grid, nb_panel, nb_edge, k_nb)

                # For conservation: M_out_p + M_out_q = 0
                # Average the two independent estimates
                avg = 0.5 * (M_out_p - M_out_q)

                # Set the matched fluxes
                _set_boundary_mass_flux!(flux_xi, flux_eta, grid, p, edge_p, k, avg)
                _set_boundary_mass_flux!(flux_xi, flux_eta, grid, nb_panel, nb_edge, k_nb, -avg)
            end
        end
    end
end
