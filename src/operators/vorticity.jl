"""
FV3 vorticity operators on the D-grid cubed sphere.

Two variants are provided:

1. **Corner-point vorticity** (`fv_vorticity`): Vorticity at cell corners
   (vertices), computed by interpolating the exact cell-mean values to
   corner positions (average of 4 surrounding cells). This is the natural
   staggering for the momentum equation's potential vorticity flux.
   Output shape: (6, Nc+1, Nc+1) with Corner staggering.

2. **Cell-mean vorticity** (`fv_vorticity_cellmean`): Computes area-averaged
   vorticity at cell centers via Stokes' theorem directly from D-grid winds
   (circulation around the primal cell boundary / cell area). This is exact
   (no approximation beyond cell averaging) and is the foundation for the
   corner-point variant.
   Output shape: (6, Nc, Nc) with CellCenter staggering.

D-grid wind convention (normalized covariant, matching FV3):
- `u_d[p, i, j]` at UEdge (ξ_{i-1/2}, η_j): V·ê_η (projection onto unit η-vector)
- `v_d[p, i, j]` at VEdge (ξ_i, η_{j-1/2}): V·ê_ξ (projection onto unit ξ-vector)

Edge length convention:
- `dx[p, i, j]` at UEdge(i,j): length of the ξ-constant edge (runs in η-direction)
- `dy[p, i, j]` at VEdge(i,j): length of the η-constant edge (runs in ξ-direction)

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Sections 3, 5, 6.
"""

# =========================================================================
# Corner-point vorticity (interpolated from cell-mean)
# =========================================================================

"""
    fv_vorticity!(omega, u_d, v_d, grid)

Compute corner-point relative vorticity by interpolating exact cell-mean
values to corner positions.

For each interior corner (i,j), the vorticity is the average of the four
surrounding cell-mean values:

    ζ_corner[i,j] = (ω[i-1,j-1] + ω[i,j-1] + ω[i-1,j] + ω[i,j]) / 4

where ω is the exact cell-mean vorticity from Stokes' theorem.

Output has Corner staggering: shape (6, Nc+1, Nc+1). Interior corners
(i ∈ 2:Nc, j ∈ 2:Nc) are computed. Boundary corners are set to zero;
accurate boundary values require cross-panel communication.

Note: FV3's production corner vorticity uses C-grid contravariant winds
with higher-order D→A→C interpolation. This simplified version is
second-order accurate and sufficient for most applications.
"""
function fv_vorticity!(omega, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    # Compute exact cell-mean vorticity
    omega_cell = fv_vorticity_cellmean(u_d, v_d, grid)

    # Zero out all values (boundary corners will remain zero)
    fill!(omega, 0)

    # Interpolate to corners by averaging 4 surrounding cell-mean values
    for p in 1:6, i in 2:Nc, j in 2:Nc
        omega[p, i, j] = 0.25 * (omega_cell[p, i - 1, j - 1] + omega_cell[p, i, j - 1] +
                                  omega_cell[p, i - 1, j] + omega_cell[p, i, j])
    end

    return omega
end

"""
    fv_vorticity(u_d, v_d, grid)

Allocating version of corner-point `fv_vorticity!`.
Returns omega array with Corner staggering (6, Nc+1, Nc+1).
"""
function fv_vorticity(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    omega = zeros(T, 6, Nc + 1, Nc + 1)
    fv_vorticity!(omega, u_d, v_d, grid)
    return omega
end

# =========================================================================
# Cell-mean vorticity (Stokes' theorem on primal cells)
# =========================================================================

"""
    fv_vorticity_cellmean!(omega, u_d, v_d, grid)

Compute cell-mean relative vorticity at cell centers using Stokes' theorem.

The circulation integral ∮ V·dl around the primal cell boundary uses the
D-grid tangential wind components times their physical edge lengths:

    Γ = v_d[i, j] · dy[i, j]           (south edge, +ξ direction)
      + u_d[i+1, j] · dx[i+1, j]       (east edge, +η direction)
      - v_d[i, j+1] · dy[i, j+1]       (north edge, -ξ direction)
      - u_d[i, j] · dx[i, j]           (west edge, -η direction)

Then: ω = Γ / A

This is exact: no finite-difference approximation is involved, just the
midpoint rule for edge integrals.

Output has CellCenter staggering: shape (6, Nc, Nc).
"""
function fv_vorticity_cellmean!(omega, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    for p in 1:6, i in 1:Nc, j in 1:Nc
        circ_south = v_d[p, i, j] * grid.dy[p, i, j]
        circ_east = u_d[p, i + 1, j] * grid.dx[p, i + 1, j]
        circ_north = -v_d[p, i, j + 1] * grid.dy[p, i, j + 1]
        circ_west = -u_d[p, i, j] * grid.dx[p, i, j]

        omega[p, i, j] = (circ_south + circ_east + circ_north + circ_west) / grid.area[p, i, j]
    end

    return omega
end

"""
    fv_vorticity_cellmean(u_d, v_d, grid)

Allocating version of `fv_vorticity_cellmean!`. Returns omega array (6, Nc, Nc).
"""
function fv_vorticity_cellmean(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    omega = zeros(T, 6, Nc, Nc)
    fv_vorticity_cellmean!(omega, u_d, v_d, grid)
    return omega
end

"""
    fv_vorticity_cellmean_arrayop(u_d, v_d, grid)

ArrayOp version of cell-mean vorticity. Returns an ArrayOp{SymReal} of shape [6, Nc, Nc].
"""
function fv_vorticity_cellmean_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    u_c = const_wrap(unwrap(u_d)); v_c = const_wrap(unwrap(v_d))
    A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)

    circ = wrap(v_c[p, i, j]) * wrap(dy_c[p, i, j]) +
           wrap(u_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) -
           wrap(v_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) -
           wrap(u_c[p, i, j]) * wrap(dx_c[p, i, j])

    expr = circ / wrap(A_c[p, i, j])
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

# =========================================================================
# Absolute vorticity (cell-mean + Coriolis)
# =========================================================================

"""
    fv_absolute_vorticity!(omega_abs, u_d, v_d, grid; Omega_rot=7.292e-5)

Compute cell-mean absolute vorticity: ω_abs = ω_rel + f

where f = 2Ω sin(lat) is the Coriolis parameter and Ω is the planetary
rotation rate. Uses cell-mean vorticity (CellCenter staggering).
"""
function fv_absolute_vorticity!(omega_abs, u_d, v_d, grid::CubedSphereGrid;
                                 Omega_rot = 7.292e-5)
    fv_vorticity_cellmean!(omega_abs, u_d, v_d, grid)
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
