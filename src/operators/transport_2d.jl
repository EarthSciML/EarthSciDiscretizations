"""
2D transport operators on the cubed-sphere grid.

Provides both the original Lax-Friedrichs ArrayOp operator (for symbolic use)
and a Lin-Rood dimensionally-split PPM operator (for accurate transport).
"""

"""
    transport_2d(q, courant_xi, courant_eta, grid)

ArrayOp for 2D Lax-Friedrichs transport tendency [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).

Note: This is a first-order scheme retained for backward compatibility.
For accurate transport, use `transport_2d_linrood!`.
"""
function transport_2d(q, courant_xi, courant_eta, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    q_c = const_wrap(unwrap(q)); cxi = const_wrap(unwrap(courant_xi)); ceta = const_wrap(unwrap(courant_eta))
    A_c = const_wrap(grid.area); dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)

    # Index (p,i,j) maps to physical cell (p, i+1, j+1)
    F_e_a = wrap(cxi[p, i + 2, j + 1]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 2, j + 1])) / 2
    F_e_d = abs(wrap(cxi[p, i + 2, j + 1])) * (wrap(q_c[p, i + 2, j + 1]) - wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_a = wrap(cxi[p, i + 1, j + 1]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_d = abs(wrap(cxi[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i, j + 1])) / 2
    G_n_a = wrap(ceta[p, i + 1, j + 2]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 1, j + 2])) / 2
    G_n_d = abs(wrap(ceta[p, i + 1, j + 2])) * (wrap(q_c[p, i + 1, j + 2]) - wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_a = wrap(ceta[p, i + 1, j + 1]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_d = abs(wrap(ceta[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i + 1, j])) / 2

    F_east = F_e_a - F_e_d; F_west = F_w_a - F_w_d
    G_north = G_n_a - G_n_d; G_south = G_s_a - G_s_d
    # Use proper FV divergence with physical edge lengths and cell areas
    xi_contrib = -(F_east * wrap(dx_c[p, i + 2, j + 1]) - F_west * wrap(dx_c[p, i + 1, j + 1])) / wrap(A_c[p, i + 1, j + 1])
    eta_contrib = -(G_north * wrap(dy_c[p, i + 1, j + 2]) - G_south * wrap(dy_c[p, i + 1, j + 1])) / wrap(A_c[p, i + 1, j + 1])
    expr = xi_contrib + eta_contrib
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end

"""
    _advective_tendency!(tend_adv, tend_flux, q, vel, grid, dim)

Convert a flux-form tendency to advective-form by adding the convergence
deformation correction:
    (dq/dt)_advective = (dq/dt)_flux + q · C_def

where C_def = (1/A) · δ(v · edge_length) is the velocity convergence.
The advective form correctly preserves uniform mixing ratios in
convergent/divergent flow fields.
"""
function _advective_tendency!(tend_adv, tend_flux, q, vel, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        for p in 1:6, i in 1:Nc, j in 1:Nc
            c_def = (vel[p, i + 1, j] * grid.dx[p, i + 1, j] -
                     vel[p, i, j] * grid.dx[p, i, j]) / grid.area[p, i, j]
            tend_adv[p, i, j] = tend_flux[p, i, j] + q[p, i, j] * c_def
        end
    else  # :eta
        for p in 1:6, i in 1:Nc, j in 1:Nc
            c_def = (vel[p, i, j + 1] * grid.dy[p, i, j + 1] -
                     vel[p, i, j] * grid.dy[p, i, j]) / grid.area[p, i, j]
            tend_adv[p, i, j] = tend_flux[p, i, j] + q[p, i, j] * c_def
        end
    end
    return tend_adv
end

"""
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)

Lin-Rood (1996) dimensionally-split 2D transport with PPM reconstruction.

Implements the dimensionally-split scheme:
    q^{n+1} = q^n + F[u*, dt, q^θ] + G[v*, dt, q^λ]

where the inner operators use the ADVECTIVE form for the half-step updates
(preserving uniform mixing ratios) and the FLUX form for the full step
(ensuring conservation):
    q^θ = q^n + (1/2)·g[v*, dt, q^n]     (advective half-step in η)
    q^λ = q^n + (1/2)·f[u*, dt, q^n]     (advective half-step in ξ)

This eliminates first-order splitting errors present in naive operator splitting.

Arguments:
- `tendency`: output array (6, Nc, Nc), modified in-place
- `q`: scalar field (6, Nc, Nc)
- `vel_xi`: ξ-velocity at cell edges (6, Nc+1, Nc)
- `vel_eta`: η-velocity at cell edges (6, Nc, Nc+1)
- `grid`: CubedSphereGrid
- `dt`: time step
"""
function transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid::CubedSphereGrid, dt)
    Nc = grid.Nc

    # Temporary arrays
    tend_xi = zeros(6, Nc, Nc)
    tend_eta = zeros(6, Nc, Nc)
    tend_adv = zeros(6, Nc, Nc)
    q_theta = zeros(6, Nc, Nc)
    q_lambda = zeros(6, Nc, Nc)

    # Step 1: Half-step intermediate in η-direction using ADVECTIVE form
    # q^θ = q^n + (1/2)·g[v*, dt, q^n]
    flux_1d_ppm!(tend_eta, q, vel_eta, grid, :eta, dt)
    _advective_tendency!(tend_adv, tend_eta, q, vel_eta, grid, :eta)
    @. q_theta = q + 0.5 * dt * tend_adv

    # Step 2: Half-step intermediate in ξ-direction using ADVECTIVE form
    # q^λ = q^n + (1/2)·f[u*, dt, q^n]
    flux_1d_ppm!(tend_xi, q, vel_xi, grid, :xi, dt)
    _advective_tendency!(tend_adv, tend_xi, q, vel_xi, grid, :xi)
    @. q_lambda = q + 0.5 * dt * tend_adv

    # Step 3: Full step with FLUX form (for conservation)
    # Compute raw fluxes for both directions
    raw_flux_xi = zeros(6, Nc + 1, Nc)
    raw_flux_eta = zeros(6, Nc, Nc + 1)
    _compute_ppm_fluxes!(raw_flux_xi, q_theta, vel_xi, grid, :xi, dt)
    _compute_ppm_fluxes!(raw_flux_eta, q_lambda, vel_eta, grid, :eta, dt)

    # Match fluxes at rotated panel boundaries for conservation
    _match_rotated_boundary_fluxes!(raw_flux_xi, raw_flux_eta, grid)

    # Convert matched fluxes to tendencies
    _flux_to_tendency!(tend_xi, raw_flux_xi, grid, :xi)
    _flux_to_tendency!(tend_eta, raw_flux_eta, grid, :eta)

    @. tendency = tend_xi + tend_eta

    return tendency
end
