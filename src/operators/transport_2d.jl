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
    Nc = grid.Nc; dξ = grid.dξ; dη = grid.dη
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    q_c = const_wrap(unwrap(q)); cxi = const_wrap(unwrap(courant_xi)); ceta = const_wrap(unwrap(courant_eta))

    F_e_a = wrap(cxi[p, i + 2, j + 1]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 2, j + 1])) / 2
    F_e_d = abs(wrap(cxi[p, i + 2, j + 1])) * (wrap(q_c[p, i + 2, j + 1]) - wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_a = wrap(cxi[p, i + 1, j + 1]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_d = abs(wrap(cxi[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i, j + 1])) / 2
    G_n_a = wrap(ceta[p, i + 1, j + 2]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 1, j + 2])) / 2
    G_n_d = abs(wrap(ceta[p, i + 1, j + 2])) * (wrap(q_c[p, i + 1, j + 2]) - wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_a = wrap(ceta[p, i + 1, j + 1]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_d = abs(wrap(ceta[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i + 1, j])) / 2

    expr = -((F_e_a - F_e_d) - (F_w_a - F_w_d)) / dξ - ((G_n_a - G_n_d) - (G_s_a - G_s_d)) / dη
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end

"""
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)

Lin-Rood (1996) dimensionally-split 2D transport with PPM reconstruction.

Implements the dimensionally-split scheme:
    q^{n+1} = q^n + F[u*, dt, q^θ] + G[v*, dt, q^λ]

where the inner operators provide cross-directional half-step updates:
    q^θ = q^n + (1/2)·G[v*, dt, q^n]     (half-step in η)
    q^λ = q^n + (1/2)·F[u*, dt, q^n]     (half-step in ξ)

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
    q_theta = zeros(6, Nc, Nc)
    q_lambda = zeros(6, Nc, Nc)

    # Step 1: Half-step intermediate in η-direction
    # q^θ = q^n + (1/2)·G[v*, dt, q^n]
    flux_1d_ppm!(tend_eta, q, vel_eta, grid, :eta, dt)
    @. q_theta = q + 0.5 * dt * tend_eta

    # Step 2: Half-step intermediate in ξ-direction
    # q^λ = q^n + (1/2)·F[u*, dt, q^n]
    flux_1d_ppm!(tend_xi, q, vel_xi, grid, :xi, dt)
    @. q_lambda = q + 0.5 * dt * tend_xi

    # Step 3: Full step with cross-terms
    # tendency = F[u*, dt, q^θ] + G[v*, dt, q^λ]
    flux_1d_ppm!(tend_xi, q_theta, vel_xi, grid, :xi, dt)
    flux_1d_ppm!(tend_eta, q_lambda, vel_eta, grid, :eta, dt)

    @. tendency = tend_xi + tend_eta

    return tendency
end
