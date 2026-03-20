"""
FV3 wind operations: covariant/contravariant decomposition and C-D grid transformation.

On the non-orthogonal gnomonic cubed-sphere grid, wind vectors must be decomposed
into covariant and contravariant components. FV3 uses D-grid staggering where
prognostic winds are normalized covariant (face-tangential velocity projections
onto unit coordinate vectors), and the C-D grid algorithm converts them to
contravariant (face-normal) winds for advection.

## Wind convention (normalized covariant, matching FV3)

All D-grid winds in this module are **normalized covariant** components:
projections of the velocity vector onto local unit tangent vectors ê_ξ, ê_η.
They have units of m/s (physical velocity). This matches FV3's convention
(Harris et al. 2021, Section 3.2, Table 5.1).

## Staggering conventions

**D-grid** (prognostic, normalized covariant):
- `u_d[p, i, j]` at UEdge position (ξ_{i-1/2}, η_j): V·ê_η
  (projection onto unit η-vector at the ξ-constant edge)
- `v_d[p, i, j]` at VEdge position (ξ_i, η_{j-1/2}): V·ê_ξ
  (projection onto unit ξ-vector at the η-constant edge)

**C-grid** (advecting, contravariant):
- `uc[p, i, j]` at UEdge position: contravariant ξ-component ũ
  (face-normal velocity through the ξ-edge)
- `vc[p, i, j]` at VEdge position: contravariant η-component ṽ
  (face-normal velocity through the η-edge)

## Conversion formulas (Harris et al. 2021, Eq. 3.4-3.7)

Given normalized covariant (u, v) and contravariant (ũ, ṽ):
    ũ = (u - v·cos α) / sin²α
    ṽ = (v - u·cos α) / sin²α

where α is the angle between ê_ξ and ê_η.

The face-normal velocity through a ξ-edge (for flux computation) is:
    U_n = ũ · sin α = (u - v·cos α) / sin α

References:
- Lin & Rood (1997), MWR
- Harris et al. (2021), GFDL FV3 Technical Memorandum, Sections 3 and 6.
"""

"""
    covariant_to_contravariant(u_cov, v_cov, sin_alpha, cos_alpha)

Convert normalized covariant wind components to contravariant at a single point.

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

Convert contravariant wind components to normalized covariant at a single point.

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

Convert D-grid normalized covariant winds to C-grid contravariant winds
using the FV3 D→A→C algorithm.

**Step 1 (D→A)**: Average D-grid covariant winds to cell centers (A-grid)
using simple 2-point averaging.

**Step 2 (Contravariant conversion)**: Convert covariant to contravariant
at cell centers using the non-orthogonal formula (Eq. 3.4).

**Step 3 (Ghost extension)**: Extend A-grid contravariant winds with ghost
cells using `extend_with_ghosts_vector`, which correctly applies rotation
matrices at panel boundaries.

**Step 4 (A→C)**: Interpolate A-grid contravariant winds to C-grid face
positions using 2-point averaging from the ghost-extended arrays.

D-grid winds:
- `u_d[p, i, j]` at UEdge (ξ_{i-1/2}, η_j): normalized covariant V·ê_η
- `v_d[p, i, j]` at VEdge (ξ_i, η_{j-1/2}): normalized covariant V·ê_ξ

C-grid output:
- `uc[p, i, j]` at UEdge: contravariant ξ-component (face-normal)
- `vc[p, i, j]` at VEdge: contravariant η-component (face-normal)

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Section 6.
"""
function dgrid_to_cgrid!(uc, vc, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    # Step 1-2: Average D-grid to A-grid and convert to contravariant
    ua = zeros(6, Nc, Nc)  # contravariant ξ at cell center
    va = zeros(6, Nc, Nc)  # contravariant η at cell center

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Average covariant winds to cell center
        u_center = 0.5 * (u_d[p, i, j] + u_d[p, i + 1, j])  # V·ê_η averaged
        v_center = 0.5 * (v_d[p, i, j] + v_d[p, i, j + 1])  # V·ê_ξ averaged

        # Convert to contravariant at cell center using sin/cos(α) at position 5
        sin_a = grid.sin_sg[p, i, j, 5]
        cos_a = grid.cos_sg[p, i, j, 5]
        sin2 = sin_a^2

        # ũ_ξ = (V·ê_ξ - V·ê_η · cos α) / sin²α
        ua[p, i, j] = (v_center - u_center * cos_a) / sin2
        # ũ_η = (V·ê_η - V·ê_ξ · cos α) / sin²α
        va[p, i, j] = (u_center - v_center * cos_a) / sin2
    end

    # Step 3: Extend A-grid with ghost cells (with vector rotation at boundaries)
    ua_ext, va_ext = extend_with_ghosts_vector(ua, va, grid)
    Ng = grid.Ng

    # Step 4: Interpolate A-grid to C-grid
    # uc at UEdge(i, j): between cell (i-1, j) and cell (i, j)
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        ie = i - 1 + Ng; je = j + Ng
        uc[p, i, j] = 0.5 * (ua_ext[p, ie, je] + ua_ext[p, ie + 1, je])
    end

    # vc at VEdge(i, j): between cell (i, j-1) and cell (i, j)
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        ie = i + Ng; je = j - 1 + Ng
        vc[p, i, j] = 0.5 * (va_ext[p, ie, je] + va_ext[p, ie, je + 1])
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
