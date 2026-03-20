"""
FV3 vorticity operator: exact computation via Stokes' theorem on the D-grid.

On the D-grid, cell-mean absolute vorticity can be computed exactly using
Stokes' theorem. The circulation around a cell boundary equals the integral
of vorticity over the cell area:

    ω_cell = (1/A) ∮ V · dl

where the line integral sums the D-grid tangential wind components times
their respective edge lengths around the cell boundary.

This is a key advantage of D-grid staggering: vorticity is computed without
any finite-difference approximation, regardless of grid resolution.

D-grid wind convention (same as wind_ops.jl):
- `u_d[p, i, j]` at UEdge (ξ_{i+1/2}, η_j): covariant η-component (tangent to ξ-edge)
- `v_d[p, i, j]` at VEdge (ξ_i, η_{j+1/2}): covariant ξ-component (tangent to η-edge)

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Section 5.
"""

"""
    fv_vorticity!(omega, u_d, v_d, grid)

Compute cell-mean relative vorticity at cell centers using Stokes' theorem.

For D-grid covariant winds (u_d = V·e_η, v_d = V·e_ξ), the circulation
integral ∮ V·dl uses the computational coordinate intervals:

    Γ = v_d[i, j] · Δξ         (south edge, +ξ direction)
      + u_d[i+1, j] · Δη       (east edge, +η direction)
      - v_d[i, j+1] · Δξ       (north edge, -ξ direction)
      - u_d[i, j] · Δη         (west edge, -η direction)

Then: ω = Γ / A

This works because dl = e_ξ dξ along ξ-edges and dl = e_η dη along η-edges,
so V·dl = (V·e_ξ) dξ = v_d dξ, and the midpoint rule gives v_d · Δξ.

Note: The sign convention follows right-hand rule with outward normal.
Positive circulation = counterclockwise when viewed from outside the sphere.
"""
function fv_vorticity!(omega, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    Δξ = grid.dξ; Δη = grid.dη

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Circulation = sum of covariant wind × coordinate interval around cell boundary
        # South edge (η = η_{j-1/2}): v_d (= V·e_ξ) in +ξ direction
        circ_south = v_d[p, i, j] * Δξ
        # East edge (ξ = ξ_{i+1/2}): u_d (= V·e_η) in +η direction
        circ_east = u_d[p, i + 1, j] * Δη
        # North edge (η = η_{j+1/2}): v_d in -ξ direction
        circ_north = -v_d[p, i, j + 1] * Δξ
        # West edge (ξ = ξ_{i-1/2}): u_d in -η direction
        circ_west = -u_d[p, i, j] * Δη

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
