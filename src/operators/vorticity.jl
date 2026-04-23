"""
FV3 vorticity operators on the D-grid cubed sphere.

Two variants are provided:

1. **Corner-point vorticity** (`fv_vorticity`): Vorticity at cell corners
   (vertices), computed by interpolating the exact cell-mean values to
   corner positions (average of 4 surrounding cells). Boundary corners
   use ghost-extended cell-mean vorticity from neighboring panels via
   `extend_with_ghosts`. This is the natural staggering for the momentum
   equation's potential vorticity flux.
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

For each corner (i,j), the vorticity is the average of the four
surrounding cell-mean values:

    ζ_corner[i,j] = (ω[i-1,j-1] + ω[i,j-1] + ω[i-1,j] + ω[i,j]) / 4

where ω is the exact cell-mean vorticity from Stokes' theorem.

Output has Corner staggering: shape (6, Nc+1, Nc+1). All corners including
boundary corners are computed. Boundary corners use ghost-extended cell-mean
vorticity from neighboring panels via `extend_with_ghosts`, which handles
cross-panel scalar ghost cell fill including cube-vertex corners where three
panels meet.
"""
function fv_vorticity!(omega, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    Ng = grid.Ng

    # Compute exact cell-mean vorticity
    omega_cell = fv_vorticity_cellmean(u_d, v_d, grid)

    # Ghost-extend cell-mean vorticity for cross-panel boundary corner values
    omega_ext = extend_with_ghosts(omega_cell, grid)

    # Interpolate to ALL corners (including boundary) by averaging 4 surrounding cell-mean values
    for p in 1:6, i in 1:(Nc + 1), j in 1:(Nc + 1)
        ie = i - 1 + Ng  # extended index for cell (i-1)
        je = j - 1 + Ng  # extended index for cell (j-1)
        omega[p, i, j] = 0.25 * (
            omega_ext[p, ie, je] + omega_ext[p, ie + 1, je] +
                omega_ext[p, ie, je + 1] + omega_ext[p, ie + 1, je + 1]
        )
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
function fv_absolute_vorticity!(
        omega_abs, u_d, v_d, grid::CubedSphereGrid;
        Omega_rot = 7.292e-5
    )
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

# =========================================================================
# ArrayOp: Corner-point vorticity (ghost-extended interpolation)
# =========================================================================

"""
    fv_vorticity_arrayop(u_d, v_d, grid)

ArrayOp version of corner-point vorticity. Computes cell-mean vorticity,
ghost-extends it via `extend_with_ghosts`, then builds an ArrayOp that
interpolates to all corners (including boundary) by averaging 4 surrounding
cell-mean values from the ghost-extended array.

Returns an ArrayOp{SymReal} of shape [6, Nc+1, Nc+1].
"""
function fv_vorticity_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc; Ng = grid.Ng

    # First compute cell-mean vorticity
    omega_cell = fv_vorticity_cellmean(u_d, v_d, grid)

    # Ghost-extend
    omega_ext = extend_with_ghosts(omega_cell, grid)

    # Build ArrayOp for corner interpolation
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    o_c = const_wrap(omega_ext)
    o = Ng  # offset; corner (i,j) maps to extended cells (i-1+Ng, j-1+Ng) through (i+Ng, j+Ng)

    expr = (
        wrap(o_c[p, i - 1 + o, j - 1 + o]) + wrap(o_c[p, i + o, j - 1 + o]) +
            wrap(o_c[p, i - 1 + o, j + o]) + wrap(o_c[p, i + o, j + o])
    ) / 4

    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc + 1), j => 1:1:(Nc + 1)))
end

# =========================================================================
# ArrayOp: Absolute vorticity (cell-mean + Coriolis)
# =========================================================================

"""
    fv_absolute_vorticity_arrayop(u_d, v_d, grid; Omega_rot=7.292e-5)

ArrayOp version of cell-mean absolute vorticity: ω_abs = ω_rel + f,
where f = 2Ω sin(lat) is the Coriolis parameter.

Returns an ArrayOp{SymReal} of shape [6, Nc, Nc].
"""
function fv_absolute_vorticity_arrayop(u_d, v_d, grid::CubedSphereGrid; Omega_rot = 7.292e-5)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]

    # Cell-mean vorticity ArrayOp
    u_c = const_wrap(unwrap(u_d)); v_c = const_wrap(unwrap(v_d))
    A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)

    circ = wrap(v_c[p, i, j]) * wrap(dy_c[p, i, j]) +
        wrap(u_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) -
        wrap(v_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) -
        wrap(u_c[p, i, j]) * wrap(dx_c[p, i, j])

    omega_rel = circ / wrap(A_c[p, i, j])

    # Precompute Coriolis parameter
    f_arr = zeros(6, Nc, Nc)
    for pp in 1:6, ii in 1:Nc, jj in 1:Nc
        f_arr[pp, ii, jj] = 2 * Omega_rot * sin(grid.lat[pp, ii, jj])
    end
    f_c = const_wrap(f_arr)

    expr = omega_rel + wrap(f_c[p, i, j])
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end
