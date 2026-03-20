"""
FV3 vorticity operator: exact computation via Stokes' theorem on the D-grid.

On the D-grid, cell-mean absolute vorticity can be computed exactly using
Stokes' theorem. The circulation around a cell boundary equals the integral
of vorticity over the cell area:

    ω_cell = (1/A) ∮ V · dl

where the line integral sums the D-grid tangential wind components times
their respective physical edge lengths around the cell boundary.

This is a key advantage of D-grid staggering: vorticity is computed without
any finite-difference approximation, regardless of grid resolution.

D-grid wind convention (normalized covariant, matching FV3):
- `u_d[p, i, j]` at UEdge (ξ_{i-1/2}, η_j): V·ê_η (projection onto unit η-vector)
- `v_d[p, i, j]` at VEdge (ξ_i, η_{j-1/2}): V·ê_ξ (projection onto unit ξ-vector)

The circulation is computed as u_d * dx (edge length) rather than u_d * Δη
(coordinate interval), because the D-grid winds are normalized covariant
components (physical velocity projections, units m/s) per Harris et al. (2021).

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Section 5.
"""

"""
    fv_vorticity!(omega, u_d, v_d, grid)

Compute cell-mean relative vorticity at cell centers using Stokes' theorem.

For D-grid normalized covariant winds (u_d = V·ê_η, v_d = V·ê_ξ), the
circulation integral ∮ V·dl uses the physical edge lengths:

    Γ = v_d[i, j] · dy[i, j]           (south edge, +ξ direction)
      + u_d[i+1, j] · dx[i+1, j]       (east edge, +η direction)
      - v_d[i, j+1] · dy[i, j+1]       (north edge, -ξ direction)
      - u_d[i, j] · dx[i, j]           (west edge, -η direction)

Then: ω = Γ / A

where dx[p,i,j] is the physical length of the ξ-constant edge at UEdge(i,j)
and dy[p,i,j] is the physical length of the η-constant edge at VEdge(i,j).

Note: The sign convention follows right-hand rule with outward normal.
Positive circulation = counterclockwise when viewed from outside the sphere.
"""
function fv_vorticity!(omega, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Circulation = normalized covariant wind × physical edge length
        # South edge (η = η_{j-1/2}): v_d (= V·ê_ξ) in +ξ direction
        circ_south = v_d[p, i, j] * grid.dy[p, i, j]
        # East edge (ξ = ξ_{i+1/2}): u_d (= V·ê_η) in +η direction
        circ_east = u_d[p, i + 1, j] * grid.dx[p, i + 1, j]
        # North edge (η = η_{j+1/2}): v_d in -ξ direction
        circ_north = -v_d[p, i, j + 1] * grid.dy[p, i, j + 1]
        # West edge (ξ = ξ_{i-1/2}): u_d in -η direction
        circ_west = -u_d[p, i, j] * grid.dx[p, i, j]

        omega[p, i, j] = (circ_south + circ_east + circ_north + circ_west) / grid.area[p, i, j]
    end

    return omega
end

"""
    fv_vorticity(u_d, v_d, grid)

Allocating version of `fv_vorticity!`. Returns omega array (6, Nc, Nc).
"""
function fv_vorticity(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    omega = zeros(T, 6, Nc, Nc)
    fv_vorticity!(omega, u_d, v_d, grid)
    return omega
end

"""
    fv_vorticity_arrayop(u_d, v_d, grid)

ArrayOp version of the FV3 vorticity operator for the symbolic discretization
pipeline. Returns an ArrayOp{SymReal} of shape [6, Nc, Nc].
"""
function fv_vorticity_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    u_c = const_wrap(unwrap(u_d)); v_c = const_wrap(unwrap(v_d))
    A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)

    # Circulation: v_d * dy (south) + u_d * dx (east) - v_d * dy (north) - u_d * dx (west)
    circ = wrap(v_c[p, i, j]) * wrap(dy_c[p, i, j]) +
           wrap(u_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) -
           wrap(v_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) -
           wrap(u_c[p, i, j]) * wrap(dx_c[p, i, j])

    expr = circ / wrap(A_c[p, i, j])
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

"""
    fv_absolute_vorticity!(omega_abs, u_d, v_d, grid; Omega_rot=7.292e-5)

Compute cell-mean absolute vorticity: ω_abs = ω_rel + f

where f = 2Ω sin(lat) is the Coriolis parameter and Ω is the planetary
rotation rate.
"""
function fv_absolute_vorticity!(omega_abs, u_d, v_d, grid::CubedSphereGrid;
                                 Omega_rot = 7.292e-5)
    fv_vorticity!(omega_abs, u_d, v_d, grid)
    Nc = grid.Nc

    for p in 1:6, i in 1:Nc, j in 1:Nc
        f = 2 * Omega_rot * sin(grid.lat[p, i, j])
        omega_abs[p, i, j] += f
    end

    return omega_abs
end

"""
    fv_absolute_vorticity(u_d, v_d, grid; Omega_rot=7.292e-5)

Allocating version of `fv_absolute_vorticity!`. Returns omega_abs (6, Nc, Nc).
"""
function fv_absolute_vorticity(u_d, v_d, grid::CubedSphereGrid; Omega_rot = 7.292e-5)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    omega_abs = zeros(T, 6, Nc, Nc)
    fv_absolute_vorticity!(omega_abs, u_d, v_d, grid; Omega_rot)
    return omega_abs
end
