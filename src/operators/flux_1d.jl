"""
1D flux-form transport operators.

The Lax-Friedrichs face-flux core is mirrored declaratively as
`discretizations/finite_volume/lax_friedrichs_flux.json`. The PPM
face-flux composition is declared as
`discretizations/finite_volume/flux_1d_ppm.json` (dsc-r1i). Imperative
PPM helpers (_compute_ppm_fluxes!, _flux_to_tendency!, boundary-flux matchers,
_build_ppm_face_expr_xi/eta) retired in esd-vx3.
"""

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


# ============================================================================
# ArrayOp-based transport operators
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
