"""
FV3 wind operations: covariant/contravariant decomposition and C-D grid transformation.

On the non-orthogonal gnomonic cubed-sphere grid, wind vectors must be decomposed
into covariant and contravariant components. FV3 uses D-grid staggering where
prognostic winds are covariant (face-tangential), and the C-D grid algorithm
converts them to contravariant (face-normal) winds for advection.

## Staggering conventions

**D-grid** (prognostic, covariant):
- `u_d[p, i, j]` at UEdge position (ξ_{i+1/2}, η_j): covariant η-component v·e_η
  (tangential to the ξ-edge, which runs in the η-direction)
- `v_d[p, i, j]` at VEdge position (ξ_i, η_{j+1/2}): covariant ξ-component v·e_ξ
  (tangential to the η-edge, which runs in the ξ-direction)

**C-grid** (advecting, contravariant):
- `uc[p, i, j]` at UEdge position: contravariant ξ-component ũ
  (face-normal velocity through the ξ-edge)
- `vc[p, i, j]` at VEdge position: contravariant η-component ṽ
  (face-normal velocity through the η-edge)

## Conversion formulas (Harris et al. 2021, Eq. 3.4-3.7)

Given covariant (u, v) and contravariant (ũ, ṽ):
    ũ = (u - v·cos α) / sin²α
    ṽ = (v - u·cos α) / sin²α

where α is the angle between e_ξ and e_η.

The face-normal velocity through a ξ-edge (for flux computation) is:
    U_n = ũ · sin α = (u - v·cos α) / sin α

References:
- Lin & Rood (1997), MWR
- Harris et al. (2021), GFDL FV3 Technical Memorandum, Sections 3 and 6.
"""

"""
    covariant_to_contravariant!(u_contra, v_contra, u_cov, v_cov, sin_alpha, cos_alpha)

Convert covariant wind components to contravariant at a single point.

    ũ = (u - v·cos α) / sin²α
    ṽ = (v - u·cos α) / sin²α
"""
function covariant_to_contravariant(u_cov, v_cov, sin_alpha, cos_alpha)
    sin2 = sin_alpha^2
    u_contra = (u_cov - v_cov * cos_alpha) / sin2
    v_contra = (v_cov - u_cov * cos_alpha) / sin2
    return (u_contra, v_contra)
end

"""
    contravariant_to_covariant(u_contra, v_contra, sin_alpha, cos_alpha)

Convert contravariant wind components to covariant at a single point.

    u = ũ + ṽ·cos α
    v = ṽ + ũ·cos α
"""
function contravariant_to_covariant(u_contra, v_contra, sin_alpha, cos_alpha)
    u_cov = u_contra + v_contra * cos_alpha
    v_cov = v_contra + u_contra * cos_alpha
    return (u_cov, v_cov)
end

"""
    dgrid_to_cgrid!(uc, vc, u_d, v_d, grid)

Interpolate D-grid covariant winds to C-grid positions and convert to
contravariant (face-normal) components.

D-grid winds:
- `u_d[p, i, j]` at UEdge (ξ_{i+1/2}, η_j): stores the covariant η-component
- `v_d[p, i, j]` at VEdge (ξ_i, η_{j+1/2}): stores the covariant ξ-component

C-grid output:
- `uc[p, i, j]` at UEdge: contravariant ξ-component (face-normal)
- `vc[p, i, j]` at VEdge: contravariant η-component (face-normal)

The conversion follows FV3 (Harris et al. 2021, Eq. 6.20):
At each UEdge position (i+1/2, j), the D-grid covariant ξ-component (v_d)
is averaged from the 4 neighboring VEdge positions, and then the pair is
converted to contravariant using the local sin/cos(α).

Similarly for VEdge positions.
"""
function dgrid_to_cgrid!(uc, vc, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    # Extend D-grid winds with ghost cells for boundary averaging
    # u_d is at UEdge (Nc+1, Nc), v_d is at VEdge (Nc, Nc+1)
    # For averaging, we need ghost cells of the cross-component
    v_d_ext = extend_with_ghosts(v_d, grid, VEdge)
    u_d_ext = extend_with_ghosts(u_d, grid, UEdge)
    Ng = grid.Ng

    # At each UEdge position (i+1/2, j), compute C-grid ξ-velocity:
    # Need covariant ξ-component (v_d) averaged to UEdge, and local u_d
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        # u_d at this position is the covariant η-component
        u_cov_eta = u_d[p, i, j]

        # Average v_d (covariant ξ-component) from 4 neighboring VEdge positions
        # VEdge neighbors of UEdge(i+1/2, j): V(i-1, j-1/2), V(i, j-1/2), V(i-1, j+1/2), V(i, j+1/2)
        # In extended array: V(i-1+Ng, j-1+Ng) etc.
        ie = i - 1 + Ng; je = j + Ng  # map UEdge i to left cell index in ext array
        v_avg = 0.25 * (v_d_ext[p, ie, je] + v_d_ext[p, ie + 1, je] +
                         v_d_ext[p, ie, je + 1] + v_d_ext[p, ie + 1, je + 1])

        # Get sin/cos(α) at this UEdge position
        # UEdge(i, j) is at ξ_{i-1/2}, η_j which is the west mid-edge of cell i
        # → sin_sg position 1 (west) of cell i, or position 3 (east) of cell i-1
        if i >= 1 && i <= Nc
            sin_a = grid.sin_sg[p, i, j, 1]
            cos_a = grid.cos_sg[p, i, j, 1]
        elseif i == Nc + 1
            sin_a = grid.sin_sg[p, Nc, j, 3]
            cos_a = grid.cos_sg[p, Nc, j, 3]
        else
            sin_a = grid.sin_sg[p, 1, j, 1]
            cos_a = grid.cos_sg[p, 1, j, 1]
        end

        # Convert: covariant (v_avg=ξ-comp, u_cov_eta=η-comp) → contravariant ξ-comp
        # ũ = (u_ξ_cov - u_η_cov·cos α) / sin²α
        # Here u_ξ_cov = v_avg, u_η_cov = u_cov_eta
        sin2 = sin_a^2
        uc[p, i, j] = (v_avg - u_cov_eta * cos_a) / sin2
    end

    # At each VEdge position (i, j+1/2), compute C-grid η-velocity:
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        # v_d at this position is the covariant ξ-component
        v_cov_xi = v_d[p, i, j]

        # Average u_d (covariant η-component) from 4 neighboring UEdge positions
        ie = i + Ng; je = j - 1 + Ng
        u_avg = 0.25 * (u_d_ext[p, ie, je] + u_d_ext[p, ie + 1, je] +
                         u_d_ext[p, ie, je + 1] + u_d_ext[p, ie + 1, je + 1])

        # Get sin/cos(α) at this VEdge position
        # VEdge(i, j) is at ξ_i, η_{j-1/2} which is the south mid-edge of cell i,j
        # → sin_sg position 2 (south) of cell (i, j)
        if j >= 1 && j <= Nc
            sin_a = grid.sin_sg[p, i, j, 2]
            cos_a = grid.cos_sg[p, i, j, 2]
        elseif j == Nc + 1
            sin_a = grid.sin_sg[p, i, Nc, 4]
            cos_a = grid.cos_sg[p, i, Nc, 4]
        else
            sin_a = grid.sin_sg[p, i, 1, 2]
            cos_a = grid.cos_sg[p, i, 1, 2]
        end

        # Convert: covariant (v_cov_xi=ξ-comp, u_avg=η-comp) → contravariant η-comp
        # ṽ = (u_η_cov - u_ξ_cov·cos α) / sin²α
        sin2 = sin_a^2
        vc[p, i, j] = (u_avg - v_cov_xi * cos_a) / sin2
    end

    return (uc, vc)
end

"""
    dgrid_to_cgrid(u_d, v_d, grid)

Allocating version of `dgrid_to_cgrid!`. Returns (uc, vc).
"""
function dgrid_to_cgrid(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    uc = zeros(T, 6, Nc + 1, Nc)
    vc = zeros(T, 6, Nc, Nc + 1)
    dgrid_to_cgrid!(uc, vc, u_d, v_d, grid)
    return (uc, vc)
end

"""
    compute_flux_with_sinsg!(flux_xi, flux_eta, vel_xi, vel_eta, grid, dt)

Compute metric-corrected mass fluxes using sin_sg with upwind selection.

Following FV3 (Harris et al. 2021), the flux through a ξ-edge includes
the sin(α) factor from the appropriate side of the edge depending on
the flow direction:

    flux_ξ[i,j] = vel_ξ[i,j] · dy[i,j] · dt · sin_sg_upwind

where sin_sg_upwind selects position 3 (east mid-edge) of the upwind cell
for positive flow, or position 1 (west mid-edge) of the downwind cell
for negative flow.
"""
function compute_flux_with_sinsg!(flux_xi, flux_eta, vel_xi, vel_eta,
                                   grid::CubedSphereGrid, dt)
    Nc = grid.Nc

    # ξ-direction fluxes at UEdge positions
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        if vel_xi[p, i, j] > 0
            # Upwind is to the west (cell i-1), use its east mid-edge sin_sg (pos 3)
            i_up = clamp(i - 1, 1, Nc)
            sin_a = grid.sin_sg[p, i_up, j, 3]
        else
            # Upwind is to the east (cell i), use its west mid-edge sin_sg (pos 1)
            i_dn = clamp(i, 1, Nc)
            sin_a = grid.sin_sg[p, i_dn, j, 1]
        end
        flux_xi[p, i, j] = vel_xi[p, i, j] * grid.dx[p, i, j] * dt * sin_a
    end

    # η-direction fluxes at VEdge positions
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        if vel_eta[p, i, j] > 0
            # Upwind is to the south (cell j-1), use its north mid-edge sin_sg (pos 4)
            j_up = clamp(j - 1, 1, Nc)
            sin_a = grid.sin_sg[p, i, j_up, 4]
        else
            # Upwind is to the north (cell j), use its south mid-edge sin_sg (pos 2)
            j_dn = clamp(j, 1, Nc)
            sin_a = grid.sin_sg[p, i, j_dn, 2]
        end
        flux_eta[p, i, j] = vel_eta[p, i, j] * grid.dy[p, i, j] * dt * sin_a
    end

    return (flux_xi, flux_eta)
end

"""
    compute_flux_with_sinsg(vel_xi, vel_eta, grid, dt)

Allocating version of `compute_flux_with_sinsg!`. Returns (flux_xi, flux_eta).
"""
function compute_flux_with_sinsg(vel_xi, vel_eta, grid::CubedSphereGrid, dt)
    Nc = grid.Nc
    T = promote_type(eltype(vel_xi), eltype(vel_eta))
    flux_xi = zeros(T, 6, Nc + 1, Nc)
    flux_eta = zeros(T, 6, Nc, Nc + 1)
    compute_flux_with_sinsg!(flux_xi, flux_eta, vel_xi, vel_eta, grid, dt)
    return (flux_xi, flux_eta)
end
