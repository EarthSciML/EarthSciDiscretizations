"""
FV3 wind operations: covariant/contravariant decomposition, C-D grid transformation,
and metric-corrected flux computation with sub-grid sin(α) upwind selection.

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

## Sub-grid sin(α) positions

FV3 stores sin(α) at 5 sub-grid positions per cell (i, j):
- Position 1: west mid-edge   (ξ_{i-1/2}, η_j)
- Position 2: south mid-edge  (ξ_i, η_{j-1/2})
- Position 3: east mid-edge   (ξ_{i+1/2}, η_j)
- Position 4: north mid-edge  (ξ_i, η_{j+1/2})
- Position 5: cell center     (ξ_i, η_j)

For flux computation, the upwind sin(α) is selected from the appropriate
sub-grid position of the upwind cell. At panel boundaries where the upwind
cell is on a neighbor panel, the correct sub-grid position is the one
coincident with the interface (e.g., position 1 at the west boundary,
position 3 at the east boundary).

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

Compute metric-corrected mass fluxes using sub-grid sin(α) with upwind selection.

Following FV3 (Harris et al. 2021), the flux through a cell edge includes
the sin(α) factor evaluated at the sub-grid position on the upwind side of
that edge:

    flux_ξ[i,j] = vel_ξ[i,j] · dx[i,j] · dt · sin_sg_upwind

For interior interfaces, sin_sg_upwind selects position 3 (east mid-edge)
of cell (i-1) for positive flow, or position 1 (west mid-edge) of cell i
for negative flow.

At panel boundaries the upwind cell lies on a neighbor panel, so the
correct sin(α) is the sub-grid value at the interface itself:
- West boundary (i=1, positive flow): position 1 of cell 1
- East boundary (i=Nc+1, negative flow): position 3 of cell Nc
- South boundary (j=1, positive flow): position 2 of cell 1
- North boundary (j=Nc+1, negative flow): position 4 of cell Nc
"""
function compute_flux_with_sinsg!(flux_xi, flux_eta, vel_xi, vel_eta,
                                   grid::CubedSphereGrid, dt)
    Nc = grid.Nc

    # ξ-direction fluxes at UEdge positions
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        if vel_xi[p, i, j] > 0
            if i - 1 >= 1
                # Interior: upwind cell i-1, use its east mid-edge (pos 3)
                sin_a = grid.sin_sg[p, i - 1, j, 3]
            else
                # West boundary: interface at ξ_{1/2}, use west mid-edge of cell 1 (pos 1)
                sin_a = grid.sin_sg[p, 1, j, 1]
            end
        else
            if i <= Nc
                # Interior: upwind cell i, use its west mid-edge (pos 1)
                sin_a = grid.sin_sg[p, i, j, 1]
            else
                # East boundary: interface at ξ_{Nc+1/2}, use east mid-edge of cell Nc (pos 3)
                sin_a = grid.sin_sg[p, Nc, j, 3]
            end
        end
        flux_xi[p, i, j] = vel_xi[p, i, j] * grid.dx[p, i, j] * dt * sin_a
    end

    # η-direction fluxes at VEdge positions
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        if vel_eta[p, i, j] > 0
            if j - 1 >= 1
                # Interior: upwind cell j-1, use its north mid-edge (pos 4)
                sin_a = grid.sin_sg[p, i, j - 1, 4]
            else
                # South boundary: interface at η_{1/2}, use south mid-edge of cell 1 (pos 2)
                sin_a = grid.sin_sg[p, i, 1, 2]
            end
        else
            if j <= Nc
                # Interior: upwind cell j, use its south mid-edge (pos 2)
                sin_a = grid.sin_sg[p, i, j, 2]
            else
                # North boundary: interface at η_{Nc+1/2}, use north mid-edge of cell Nc (pos 4)
                sin_a = grid.sin_sg[p, i, Nc, 4]
            end
        end
        flux_eta[p, i, j] = vel_eta[p, i, j] * grid.dy[p, i, j] * dt * sin_a
    end

    return (flux_xi, flux_eta)
end

"""
    compute_flux_with_sinsg(vel_xi, vel_eta, grid, dt)

Allocating version of `compute_flux_with_sinsg!`. Returns `(flux_xi, flux_eta)`.
"""
function compute_flux_with_sinsg(vel_xi, vel_eta, grid::CubedSphereGrid, dt)
    Nc = grid.Nc
    T = promote_type(eltype(vel_xi), eltype(vel_eta))
    flux_xi = zeros(T, 6, Nc + 1, Nc)
    flux_eta = zeros(T, 6, Nc, Nc + 1)
    compute_flux_with_sinsg!(flux_xi, flux_eta, vel_xi, vel_eta, grid, dt)
    return (flux_xi, flux_eta)
end

# --- Precomputed sin_sg arrays for ArrayOp flux computation ---

"""
    _precompute_sinsg_xi(grid) -> (sin_pos, sin_neg)

Precompute sin(α) arrays for ξ-direction flux at each UEdge interface.

Returns two arrays of size `(6, Nc+1, Nc)`:
- `sin_pos[p,i,j]`: sin(α) to use when `vel_xi[p,i,j] > 0` (positive/eastward flow)
- `sin_neg[p,i,j]`: sin(α) to use when `vel_xi[p,i,j] < 0` (negative/westward flow)

For interior interfaces, `sin_pos` uses position 3 (east mid-edge) of cell (i-1)
and `sin_neg` uses position 1 (west mid-edge) of cell i. At panel boundaries
the values are the sub-grid sin(α) coincident with the interface.
"""
function _precompute_sinsg_xi(grid::CubedSphereGrid)
    Nc = grid.Nc
    sin_pos = zeros(6, Nc + 1, Nc)
    sin_neg = zeros(6, Nc + 1, Nc)

    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        # Positive flow: upwind is cell (i-1)
        if i - 1 >= 1
            sin_pos[p, i, j] = grid.sin_sg[p, i - 1, j, 3]
        else
            # West boundary: use west mid-edge of cell 1 (pos 1)
            sin_pos[p, i, j] = grid.sin_sg[p, 1, j, 1]
        end

        # Negative flow: upwind is cell i
        if i <= Nc
            sin_neg[p, i, j] = grid.sin_sg[p, i, j, 1]
        else
            # East boundary: use east mid-edge of cell Nc (pos 3)
            sin_neg[p, i, j] = grid.sin_sg[p, Nc, j, 3]
        end
    end

    return (sin_pos, sin_neg)
end

"""
    _precompute_sinsg_eta(grid) -> (sin_pos, sin_neg)

Precompute sin(α) arrays for η-direction flux at each VEdge interface.

Returns two arrays of size `(6, Nc, Nc+1)`:
- `sin_pos[p,i,j]`: sin(α) to use when `vel_eta[p,i,j] > 0` (positive/northward flow)
- `sin_neg[p,i,j]`: sin(α) to use when `vel_eta[p,i,j] < 0` (negative/southward flow)

For interior interfaces, `sin_pos` uses position 4 (north mid-edge) of cell (j-1)
and `sin_neg` uses position 2 (south mid-edge) of cell j. At panel boundaries
the values are the sub-grid sin(α) coincident with the interface.
"""
function _precompute_sinsg_eta(grid::CubedSphereGrid)
    Nc = grid.Nc
    sin_pos = zeros(6, Nc, Nc + 1)
    sin_neg = zeros(6, Nc, Nc + 1)

    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        # Positive flow: upwind is cell (j-1)
        if j - 1 >= 1
            sin_pos[p, i, j] = grid.sin_sg[p, i, j - 1, 4]
        else
            # South boundary: use south mid-edge of cell 1 (pos 2)
            sin_pos[p, i, j] = grid.sin_sg[p, i, 1, 2]
        end

        # Negative flow: upwind is cell j
        if j <= Nc
            sin_neg[p, i, j] = grid.sin_sg[p, i, j, 2]
        else
            # North boundary: use north mid-edge of cell Nc (pos 4)
            sin_neg[p, i, j] = grid.sin_sg[p, i, Nc, 4]
        end
    end

    return (sin_pos, sin_neg)
end

# --- ArrayOp versions of sin_sg flux computation ---

"""
    compute_flux_with_sinsg_xi_arrayop(vel_xi, grid, dt)

ArrayOp for ξ-direction metric-corrected flux at UEdge positions [6, Nc+1, Nc].

Uses an upwind-blended expression that avoids branching:

    flux = ((v + |v|) · sin_pos + (v - |v|) · sin_neg) / 2 · dx · dt

where `sin_pos` and `sin_neg` are precomputed sub-grid sin(α) arrays
for positive and negative flow directions respectively (see `_precompute_sinsg_xi`).
"""
function compute_flux_with_sinsg_xi_arrayop(vel_xi, grid::CubedSphereGrid, dt)
    Nc = grid.Nc
    sin_pos, sin_neg = _precompute_sinsg_xi(grid)

    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    vel_c = const_wrap(unwrap(vel_xi))
    dx_c = const_wrap(grid.dx)
    sp_c = const_wrap(sin_pos)
    sn_c = const_wrap(sin_neg)

    v = wrap(vel_c[p, i, j])
    flux_expr = ((v + abs(v)) * wrap(sp_c[p, i, j]) +
                 (v - abs(v)) * wrap(sn_c[p, i, j])) / 2 *
                wrap(dx_c[p, i, j]) * dt

    return make_arrayop(idx, unwrap(flux_expr), Dict(p => 1:1:6, i => 1:1:(Nc + 1), j => 1:1:Nc))
end

"""
    compute_flux_with_sinsg_eta_arrayop(vel_eta, grid, dt)

ArrayOp for η-direction metric-corrected flux at VEdge positions [6, Nc, Nc+1].

Uses an upwind-blended expression that avoids branching:

    flux = ((v + |v|) · sin_pos + (v - |v|) · sin_neg) / 2 · dy · dt

where `sin_pos` and `sin_neg` are precomputed sub-grid sin(α) arrays
for positive and negative flow directions respectively (see `_precompute_sinsg_eta`).
"""
function compute_flux_with_sinsg_eta_arrayop(vel_eta, grid::CubedSphereGrid, dt)
    Nc = grid.Nc
    sin_pos, sin_neg = _precompute_sinsg_eta(grid)

    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    vel_c = const_wrap(unwrap(vel_eta))
    dy_c = const_wrap(grid.dy)
    sp_c = const_wrap(sin_pos)
    sn_c = const_wrap(sin_neg)

    v = wrap(vel_c[p, i, j])
    flux_expr = ((v + abs(v)) * wrap(sp_c[p, i, j]) +
                 (v - abs(v)) * wrap(sn_c[p, i, j])) / 2 *
                wrap(dy_c[p, i, j]) * dt

    return make_arrayop(idx, unwrap(flux_expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc + 1)))
end

# --- ArrayOp version of D→C grid conversion ---

"""
    dgrid_to_cgrid_arrayop(u_d, v_d, grid) -> (uc_arrayop, vc_arrayop)

ArrayOp version of the D→C grid conversion for the symbolic discretization pipeline.

Steps 1–3 (D→A averaging, covariant→contravariant conversion, ghost extension with
vector rotation) are evaluated numerically. Step 4 (A→C interpolation) is expressed
as ArrayOps that index into the ghost-extended A-grid arrays as constants.

Returns a tuple of two ArrayOps:
- `uc_arrayop`: contravariant ξ-component at UEdge [6, Nc+1, Nc]
- `vc_arrayop`: contravariant η-component at VEdge [6, Nc, Nc+1]
"""
function dgrid_to_cgrid_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc; Ng = grid.Ng

    # Steps 1-2: D→A averaging + covariant→contravariant (numerical)
    ua = zeros(6, Nc, Nc)
    va = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        u_center = 0.5 * (u_d[p, i, j] + u_d[p, i + 1, j])
        v_center = 0.5 * (v_d[p, i, j] + v_d[p, i, j + 1])
        sin_a = grid.sin_sg[p, i, j, 5]
        cos_a = grid.cos_sg[p, i, j, 5]
        sin2 = sin_a^2
        ua[p, i, j] = (v_center - u_center * cos_a) / sin2
        va[p, i, j] = (u_center - v_center * cos_a) / sin2
    end

    # Step 3: Ghost extension with vector rotation (numerical)
    ua_ext, va_ext = extend_with_ghosts_vector(ua, va, grid)

    # Step 4: A→C interpolation as ArrayOps
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    o = Ng  # offset: real cell k maps to extended index k + Ng

    # uc at UEdge(i, j): between A-grid cells (i-1) and (i) in ξ-direction
    ua_c = const_wrap(ua_ext)
    uc_expr = (wrap(ua_c[p, i - 1 + o, j + o]) + wrap(ua_c[p, i + o, j + o])) / 2
    uc_ao = make_arrayop(idx, unwrap(uc_expr), Dict(p => 1:1:6, i => 1:1:(Nc + 1), j => 1:1:Nc))

    # vc at VEdge(i, j): between A-grid cells (j-1) and (j) in η-direction
    va_c = const_wrap(va_ext)
    vc_expr = (wrap(va_c[p, i + o, j - 1 + o]) + wrap(va_c[p, i + o, j + o])) / 2
    vc_ao = make_arrayop(idx, unwrap(vc_expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc + 1)))

    return (uc_ao, vc_ao)
end
