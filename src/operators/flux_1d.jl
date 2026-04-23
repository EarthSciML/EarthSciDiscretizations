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
    A_c = const_wrap(grid.area); dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    if dim == :xi
        # Index (p,i,j) maps to physical cell (p, i+1, j)
        # East face = interface i+2, West face = interface i+1
        F_e_a = wrap(c_c[p, i + 2, j]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 2, j])) / 2
        F_e_d = abs(wrap(c_c[p, i + 2, j])) * (wrap(q_c[p, i + 2, j]) - wrap(q_c[p, i + 1, j])) / 2
        F_w_a = wrap(c_c[p, i + 1, j]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i + 1, j])) / 2
        F_w_d = abs(wrap(c_c[p, i + 1, j])) * (wrap(q_c[p, i + 1, j]) - wrap(q_c[p, i, j])) / 2
        F_east = F_e_a - F_e_d
        F_west = F_w_a - F_w_d
        expr = -(F_east * wrap(dx_c[p, i + 2, j]) - F_west * wrap(dx_c[p, i + 1, j])) / wrap(A_c[p, i + 1, j])
        ranges = Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:Nc)
    else
        # Index (p,i,j) maps to physical cell (p, i, j+1)
        # North face = interface j+2, South face = interface j+1
        F_n_a = wrap(c_c[p, i, j + 2]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i, j + 2])) / 2
        F_n_d = abs(wrap(c_c[p, i, j + 2])) * (wrap(q_c[p, i, j + 2]) - wrap(q_c[p, i, j + 1])) / 2
        F_s_a = wrap(c_c[p, i, j + 1]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i, j + 1])) / 2
        F_s_d = abs(wrap(c_c[p, i, j + 1])) * (wrap(q_c[p, i, j + 1]) - wrap(q_c[p, i, j])) / 2
        F_north = F_n_a - F_n_d
        F_south = F_s_a - F_s_d
        expr = -(F_north * wrap(dy_c[p, i, j + 2]) - F_south * wrap(dy_c[p, i, j + 1])) / wrap(A_c[p, i, j + 1])
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
            for i in 1:(Nc + 1)
                # Cells adjacent to this interface (in extended array)
                ie = i + Ng  # extended array index for cell i (or edge i)
                je = j + Ng

                # PPM reconstruction for the upwind cell
                # We need cell values: ie-2, ie-1, ie, ie+1, ie+2 in extended array
                qim2 = q_ext[p, ie - 2, je]; qim1 = q_ext[p, ie - 1, je]
                qi = q_ext[p, ie, je];     qip1 = q_ext[p, ie + 1, je]
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
            for j in 1:(Nc + 1)
                ie = i + Ng; je = j + Ng

                qjm2 = q_ext[p, ie, je - 2]; qjm1 = q_ext[p, ie, je - 1]
                qj = q_ext[p, ie, je];     qjp1 = q_ext[p, ie, je + 1]
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
            tendency[p, i, j] = -(
                flux[p, i + 1, j] * grid.dx[p, i + 1, j] -
                    flux[p, i, j] * grid.dx[p, i, j]
            ) / grid.area[p, i, j]
        end
    else
        for p in 1:6, i in 1:Nc, j in 1:Nc
            tendency[p, i, j] = -(
                flux[p, i, j + 1] * grid.dy[p, i, j + 1] -
                    flux[p, i, j] * grid.dy[p, i, j]
            ) / grid.area[p, i, j]
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
For each shared edge where both panels have ξ-fluxes, average the mass fluxes
(flux × edge_length) to ensure exact conservation, then convert back.
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
                # Convert to mass fluxes for averaging
                # Both fluxes are positive in the +ξ direction:
                #   Panel p east: positive = outflow from p
                #   Neighbor west: positive = inflow to neighbor
                # For conservation: outflow from p = inflow to neighbor
                dx_p = grid.dx[p, Nc + 1, j]
                dx_nb = grid.dx[nb.neighbor_panel, 1, j_nb]
                M_p = flux[p, Nc + 1, j] * dx_p
                M_nb = flux[nb.neighbor_panel, 1, j_nb] * dx_nb
                avg_M = 0.5 * (M_p + M_nb)
                flux[p, Nc + 1, j] = avg_M / dx_p
                flux[nb.neighbor_panel, 1, j_nb] = avg_M / dx_nb
            end
        end
    end
    return
end

"""
    _match_boundary_fluxes_eta!(flux, grid)

Match η-direction fluxes at panel boundaries for same-direction connections.
Uses mass fluxes (flux × edge_length) for averaging to ensure exact conservation.
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
                # Convert to mass fluxes for averaging
                dy_p = grid.dy[p, i, Nc + 1]
                dy_nb = grid.dy[nb.neighbor_panel, i_nb, 1]
                M_p = flux[p, i, Nc + 1] * dy_p
                M_nb = flux[nb.neighbor_panel, i_nb, 1] * dy_nb
                avg_M = 0.5 * (M_p + M_nb)
                flux[p, i, Nc + 1] = avg_M / dy_p
                flux[nb.neighbor_panel, i_nb, 1] = avg_M / dy_nb
            end
        end
    end
    return
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
    return if edge == East
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
    processed = Set{Tuple{Int, EdgeDirection, Int, EdgeDirection}}()

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
    return
end

# ============================================================================
# ArrayOp-based PPM transport operators
# ============================================================================

"""
    compute_courant_numbers(vel, dt, grid, dim)

Precompute Courant numbers at all cell interfaces for use with PPM ArrayOps.
Returns an array of size (6, Nc+1, Nc) for `:xi` or (6, Nc, Nc+1) for `:eta`.
"""
function compute_courant_numbers(vel, dt, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        courant = zeros(6, Nc + 1, Nc)
        for p in 1:6, i in 1:(Nc + 1), j in 1:Nc
            courant[p, i, j] = _get_courant_xi(vel, dt, grid, p, i, j)
        end
    else
        courant = zeros(6, Nc, Nc + 1)
        for p in 1:6, i in 1:Nc, j in 1:(Nc + 1)
            courant[p, i, j] = _get_courant_eta(vel, dt, grid, p, i, j)
        end
    end
    return courant
end

"""
    _build_interface_distances(grid, dim)

Build a single array of center-to-center distances at all cell interfaces,
combining interior distances and cross-panel boundary distances.

Returns (6, Nc+1, Nc) for `:xi` or (6, Nc, Nc+1) for `:eta`.
"""
function _build_interface_distances(grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        dist = zeros(6, Nc + 1, Nc)
        for p in 1:6, j in 1:Nc
            dist[p, 1, j] = grid.dist_xi_bnd[p, 1, j]
            for i in 2:Nc
                dist[p, i, j] = grid.dist_xi[p, i - 1, j]
            end
            dist[p, Nc + 1, j] = grid.dist_xi_bnd[p, 2, j]
        end
        return dist
    else
        dist = zeros(6, Nc, Nc + 1)
        for p in 1:6, i in 1:Nc
            dist[p, i, 1] = grid.dist_eta_bnd[p, i, 1]
            for j in 2:Nc
                dist[p, i, j] = grid.dist_eta[p, i, j - 1]
            end
            dist[p, i, Nc + 1] = grid.dist_eta_bnd[p, i, 2]
        end
        return dist
    end
end

"""
    compute_courant_numbers_arrayop(vel, dt, grid, dim)

ArrayOp for Courant numbers at all cell interfaces [6, Nc+1, Nc] for `:xi`
or [6, Nc, Nc+1] for `:eta`.

Uses precomputed interface distances to express the Courant number as a single
ArrayOp: `courant[p,i,j] = vel[p,i,j] * dt / dist[p,i,j]`.
"""
function compute_courant_numbers_arrayop(vel, dt, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    dist = _build_interface_distances(grid, dim)
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    v_c = const_wrap(unwrap(vel))
    d_c = const_wrap(dist)
    expr = wrap(v_c[p, i, j]) * dt / wrap(d_c[p, i, j])
    if dim == :xi
        ranges = Dict(p => 1:1:6, i => 1:1:(Nc + 1), j => 1:1:Nc)
    else
        ranges = Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc + 1))
    end
    return make_arrayop(idx, unwrap(expr), ranges)
end

"""
    _build_ppm_face_expr_xi(q_c, c_c, v_c, p, i_face, j, o)

Build a symbolic expression for the PPM flux at a single ξ-interface.
Uses `ifelse` for CW84 limiting and upwind selection, making it compatible
with ArrayOp symbolic tracing.

Arguments:
- `q_c`: Const-wrapped ghost-extended scalar field
- `c_c`: Const-wrapped Courant number at interfaces (6, Nc+1, Nc)
- `v_c`: Const-wrapped velocity at interfaces (6, Nc+1, Nc)
- `p, i_face, j`: symbolic index expressions for interface position
- `o`: ghost offset (Ng)
"""
function _build_ppm_face_expr_xi(q_c, c_c, v_c, p, i_face, j, o)
    # Interface i_face sits between cells (i_face-1) and i_face in interior coords
    # In extended array: cell k -> index k+o
    ie = i_face + o   # extended index for right cell (i_face)
    je = j + o        # extended j index

    # Stencil values from extended array (7 points)
    qm3 = wrap(q_c[p, ie - 3, je])
    qm2 = wrap(q_c[p, ie - 2, je])
    qm1 = wrap(q_c[p, ie - 1, je])  # left cell (i_face - 1)
    q0 = wrap(q_c[p, ie, je])       # right cell (i_face)
    qp1 = wrap(q_c[p, ie + 1, je])
    qp2 = wrap(q_c[p, ie + 2, je])

    # 4th-order interface value between qm1 and q0
    qi_half = (7.0 / 12.0) * (qm1 + q0) - (1.0 / 12.0) * (qm2 + qp1)

    # Left cell (qm1) reconstruction: edges ql_L, qr_L
    ql_L = (7.0 / 12.0) * (qm2 + qm1) - (1.0 / 12.0) * (qm3 + q0)
    qr_L = qi_half
    ql_L, qr_L = _ppm_limit_cw84_sym(ql_L, qr_L, qm1)

    # Right cell (q0) reconstruction: edges ql_R, qr_R
    ql_R = qi_half
    qr_R = (7.0 / 12.0) * (q0 + qp1) - (1.0 / 12.0) * (qm1 + qp2)
    ql_R, qr_R = _ppm_limit_cw84_sym(ql_R, qr_R, q0)

    # Courant number and velocity at this interface
    c = wrap(c_c[p, i_face, j])
    v = wrap(v_c[p, i_face, j])
    c_abs = abs(c)

    # PPM flux integral for left cell (positive flow, c >= 0)
    dq_L = qr_L - ql_L
    q6_L = 6.0 * (qm1 - 0.5 * (ql_L + qr_L))
    int_left = qr_L - 0.5 * c_abs * (dq_L - q6_L * (1.0 - (2.0 / 3.0) * c_abs))

    # PPM flux integral for right cell (negative flow, c < 0)
    dq_R = qr_R - ql_R
    q6_R = 6.0 * (q0 - 0.5 * (ql_R + qr_R))
    int_right = ql_R + 0.5 * c_abs * (dq_R + q6_R * (1.0 - (2.0 / 3.0) * c_abs))

    # Select based on flow direction
    flux_val = ifelse(c >= 0, int_left, int_right)

    return v * flux_val
end

"""
    _build_ppm_face_expr_eta(q_c, c_c, v_c, p, i, j_face, o)

Build a symbolic expression for the PPM flux at a single η-interface.
Same as `_build_ppm_face_expr_xi` but for the η-direction.
"""
function _build_ppm_face_expr_eta(q_c, c_c, v_c, p, i, j_face, o)
    ie = i + o
    je = j_face + o   # extended index for top cell (j_face)

    # Stencil values from extended array
    qm3 = wrap(q_c[p, ie, je - 3])
    qm2 = wrap(q_c[p, ie, je - 2])
    qm1 = wrap(q_c[p, ie, je - 1])  # bottom cell (j_face - 1)
    q0 = wrap(q_c[p, ie, je])       # top cell (j_face)
    qp1 = wrap(q_c[p, ie, je + 1])
    qp2 = wrap(q_c[p, ie, je + 2])

    qi_half = (7.0 / 12.0) * (qm1 + q0) - (1.0 / 12.0) * (qm2 + qp1)

    ql_L = (7.0 / 12.0) * (qm2 + qm1) - (1.0 / 12.0) * (qm3 + q0)
    qr_L = qi_half
    ql_L, qr_L = _ppm_limit_cw84_sym(ql_L, qr_L, qm1)

    ql_R = qi_half
    qr_R = (7.0 / 12.0) * (q0 + qp1) - (1.0 / 12.0) * (qm1 + qp2)
    ql_R, qr_R = _ppm_limit_cw84_sym(ql_R, qr_R, q0)

    c = wrap(c_c[p, i, j_face])
    v = wrap(v_c[p, i, j_face])
    c_abs = abs(c)

    dq_L = qr_L - ql_L
    q6_L = 6.0 * (qm1 - 0.5 * (ql_L + qr_L))
    int_left = qr_L - 0.5 * c_abs * (dq_L - q6_L * (1.0 - (2.0 / 3.0) * c_abs))

    dq_R = qr_R - ql_R
    q6_R = 6.0 * (q0 - 0.5 * (ql_R + qr_R))
    int_right = ql_R + 0.5 * c_abs * (dq_R + q6_R * (1.0 - (2.0 / 3.0) * c_abs))

    flux_val = ifelse(c >= 0, int_left, int_right)
    return v * flux_val
end

"""
    flux_1d_ppm_arrayop(q_ext, courant, vel, grid, dim)

ArrayOp for PPM-based 1D flux-form transport tendency [6, Nc, Nc].

Combines PPM reconstruction, CW84 limiting (via `ifelse`), upwind flux selection,
and FV divergence into a single ArrayOp expression. Operates on a ghost-extended
scalar field `q_ext` with precomputed Courant numbers and velocities at interfaces.

Arguments:
- `q_ext`: ghost-extended scalar field (6, Nc+2Ng, Nc+2Ng)
- `courant`: Courant numbers at interfaces (6, Nc+1, Nc) for :xi or (6, Nc, Nc+1) for :eta
- `vel`: velocity at interfaces, same size as courant
- `grid`: CubedSphereGrid
- `dim`: :xi or :eta
"""
function flux_1d_ppm_arrayop(q_ext, courant, vel, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc; Ng = grid.Ng; o = Ng
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]

    q_c = const_wrap(q_ext)
    c_c = const_wrap(unwrap(courant))
    v_c = const_wrap(unwrap(vel))
    A_c = const_wrap(grid.area)

    if dim == :xi
        dx_c = const_wrap(grid.dx)
        # East face = interface i+1, West face = interface i
        F_east = _build_ppm_face_expr_xi(q_c, c_c, v_c, p, i + 1, j, o)
        F_west = _build_ppm_face_expr_xi(q_c, c_c, v_c, p, i, j, o)
        expr = -(F_east * wrap(dx_c[p, i + 1, j]) - F_west * wrap(dx_c[p, i, j])) / wrap(A_c[p, i, j])
    else
        dy_c = const_wrap(grid.dy)
        # North face = interface j+1, South face = interface j
        F_north = _build_ppm_face_expr_eta(q_c, c_c, v_c, p, i, j + 1, o)
        F_south = _build_ppm_face_expr_eta(q_c, c_c, v_c, p, i, j, o)
        expr = -(F_north * wrap(dy_c[p, i, j + 1]) - F_south * wrap(dy_c[p, i, j])) / wrap(A_c[p, i, j])
    end

    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

"""
    flux_to_tendency_arrayop(flux, grid, dim)

ArrayOp for converting interface fluxes to cell tendencies via FV divergence [6, Nc, Nc].

    tendency[p,i,j] = -(F_{i+1} * edge_length_{i+1} - F_i * edge_length_i) / area
"""
function flux_to_tendency_arrayop(flux, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    flux_c = const_wrap(unwrap(flux))
    A_c = const_wrap(grid.area)

    if dim == :xi
        dx_c = const_wrap(grid.dx)
        expr = -(
            wrap(flux_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) -
                wrap(flux_c[p, i, j]) * wrap(dx_c[p, i, j])
        ) / wrap(A_c[p, i, j])
    else
        dy_c = const_wrap(grid.dy)
        expr = -(
            wrap(flux_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) -
                wrap(flux_c[p, i, j]) * wrap(dy_c[p, i, j])
        ) / wrap(A_c[p, i, j])
    end

    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

"""
    advective_tendency_arrayop(tend_flux, q, vel, grid, dim)

ArrayOp for advective-form tendency [6, Nc, Nc].

Converts a flux-form tendency to advective form by adding the velocity
convergence correction: tend_adv = tend_flux + q · C_def
where C_def = δ(v · edge_length) / area.
"""
function advective_tendency_arrayop(tend_flux, q, vel, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    tf_c = const_wrap(unwrap(tend_flux))
    q_c = const_wrap(unwrap(q))
    v_c = const_wrap(unwrap(vel))
    A_c = const_wrap(grid.area)

    if dim == :xi
        dx_c = const_wrap(grid.dx)
        c_def = (
            wrap(v_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) -
                wrap(v_c[p, i, j]) * wrap(dx_c[p, i, j])
        ) / wrap(A_c[p, i, j])
    else
        dy_c = const_wrap(grid.dy)
        c_def = (
            wrap(v_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) -
                wrap(v_c[p, i, j]) * wrap(dy_c[p, i, j])
        ) / wrap(A_c[p, i, j])
    end

    expr = wrap(tf_c[p, i, j]) + wrap(q_c[p, i, j]) * c_def
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end
