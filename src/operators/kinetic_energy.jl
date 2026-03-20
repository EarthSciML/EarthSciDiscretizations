"""
FV3 kinetic energy computation on the non-orthogonal cubed-sphere grid.

## Primary formulation: upstream-biased KE (Harris et al. 2021, Eq. 6.3)

FV3 computes KE using the mixed covariant-contravariant form:

    KE = (1/2)(u·ṽ + v·ũ)

where u = V·ê_η and v = V·ê_ξ are the D-grid normalized covariant winds,
ṽ = (u - v·cos α)/sin²α is the contravariant η-component, and
ũ = (v - u·cos α)/sin²α is the contravariant ξ-component.

The products are evaluated at D-grid face positions and averaged to cell centers:

    KE[i,j] = (1/2) · [(u_d[i,j]·ṽ[i,j] + u_d[i+1,j]·ṽ[i+1,j]) / 2 +
                         (v_d[i,j]·ũ[i,j] + v_d[i,j+1]·ũ[i,j+1]) / 2]

where ṽ is interpolated from A-grid to UEdge positions (pairing with u_d),
and ũ is interpolated from A-grid to VEdge positions (pairing with v_d).

This formulation provides upstream-biasing because the A-grid contravariant
winds are averaged to face positions rather than evaluated at cell centers,
avoiding the spurious grid-scale noise that causes the Hollingsworth-Kallberg
instability.

## Alternative: cell-center KE

A simpler cell-center formulation is also provided (functions with `_cell`
suffix) that averages D-grid winds to cell centers and uses the non-orthogonal
formula with center sin/cos(α):

    KE = (u² + v² - 2uv·cos α) / (2·sin²α)

This is equivalent to the mixed form but evaluates everything at the cell
center, which can trigger the HK instability in dynamical core simulations.

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Eq. 3.5, 6.3.
"""

# =============================================================================
# Upstream-biased KE (primary formulation)
# =============================================================================

"""
    fv_kinetic_energy!(ke, u_d, v_d, grid)

Compute cell-center kinetic energy using the FV3 upstream-biased formulation
(Harris et al. 2021, Eq. 6.3).

Uses the mixed covariant-contravariant form KE = (1/2)(u·ṽ + v·ũ) with the
products evaluated at D-grid face positions rather than cell centers, which
avoids the Hollingsworth-Kallberg instability.

Steps:
1. Average D-grid covariant winds to A-grid cell centers
2. Convert to contravariant at cell centers: ũ_ξ and ṽ_η
3. Ghost-extend A-grid contravariant fields (with vector rotation)
4. Interpolate ṽ_η to UEdge positions → pair with u_d (V·ê_η)
   Interpolate ũ_ξ to VEdge positions → pair with v_d (V·ê_ξ)

    KE[i,j] = (1/2) · [(u_d[i,j]·ṽ[i,j] + u_d[i+1,j]·ṽ[i+1,j]) / 2 +
                         (v_d[i,j]·ũ[i,j] + v_d[i,j+1]·ũ[i,j+1]) / 2]

Arguments:
- `ke`: output array (6, Nc, Nc), modified in-place
- `u_d`: D-grid η-component at UEdge (6, Nc+1, Nc): V·ê_η
- `v_d`: D-grid ξ-component at VEdge (6, Nc, Nc+1): V·ê_ξ
- `grid`: CubedSphereGrid
"""
function fv_kinetic_energy!(ke, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    # Steps 1-2: Average D-grid to A-grid and convert to contravariant
    ua = zeros(6, Nc, Nc)  # ũ_ξ (contravariant ξ) at cell center
    va = zeros(6, Nc, Nc)  # ṽ_η (contravariant η) at cell center

    for p in 1:6, i in 1:Nc, j in 1:Nc
        u_center = 0.5 * (u_d[p, i, j] + u_d[p, i + 1, j])  # V·ê_η averaged
        v_center = 0.5 * (v_d[p, i, j] + v_d[p, i, j + 1])  # V·ê_ξ averaged

        sin_a = grid.sin_sg[p, i, j, 5]
        cos_a = grid.cos_sg[p, i, j, 5]
        sin2 = sin_a^2

        ua[p, i, j] = (v_center - u_center * cos_a) / sin2  # ũ_ξ
        va[p, i, j] = (u_center - v_center * cos_a) / sin2  # ṽ_η
    end

    # Step 3: Ghost-extend A-grid contravariant (with vector rotation)
    ua_ext, va_ext = extend_with_ghosts_vector(ua, va, grid)
    Ng = grid.Ng

    # Step 4: Evaluate KE = (1/2)(u·ṽ + v·ũ) at face positions
    for p in 1:6, i in 1:Nc, j in 1:Nc
        # η-contribution: u_d · ṽ_η at UEdge positions
        # UEdge(i,j) sits between A-grid cells (i-1,j) and (i,j)
        ie_w = i - 1 + Ng; je = j + Ng
        va_at_west = 0.5 * (va_ext[p, ie_w, je] + va_ext[p, ie_w + 1, je])
        ie_e = i + Ng
        va_at_east = 0.5 * (va_ext[p, ie_e, je] + va_ext[p, ie_e + 1, je])
        ke_eta = 0.5 * (u_d[p, i, j] * va_at_west + u_d[p, i + 1, j] * va_at_east)

        # ξ-contribution: v_d · ũ_ξ at VEdge positions
        # VEdge(i,j) sits between A-grid cells (i,j-1) and (i,j)
        ie = i + Ng; je_s = j - 1 + Ng
        ua_at_south = 0.5 * (ua_ext[p, ie, je_s] + ua_ext[p, ie, je_s + 1])
        je_n = j + Ng
        ua_at_north = 0.5 * (ua_ext[p, ie, je_n] + ua_ext[p, ie, je_n + 1])
        ke_xi = 0.5 * (v_d[p, i, j] * ua_at_south + v_d[p, i, j + 1] * ua_at_north)

        ke[p, i, j] = 0.5 * (ke_xi + ke_eta)
    end

    return ke
end

"""
    fv_kinetic_energy(u_d, v_d, grid)

Allocating version of the upstream-biased `fv_kinetic_energy!`.
Returns ke array (6, Nc, Nc).
"""
function fv_kinetic_energy(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    ke = zeros(T, 6, Nc, Nc)
    fv_kinetic_energy!(ke, u_d, v_d, grid)
    return ke
end

"""
    fv_kinetic_energy_arrayop(u_d, v_d, grid)

ArrayOp version of the FV3 upstream-biased kinetic energy operator for the
symbolic discretization pipeline. Returns an ArrayOp{SymReal} of shape
[6, Nc, Nc].

Steps 1-3 (D→A averaging, covariant→contravariant conversion, ghost extension)
are evaluated numerically. Step 4 (face-position KE evaluation) is expressed
as an ArrayOp that indexes into the ghost-extended A-grid arrays and D-grid
arrays as constants.
"""
function fv_kinetic_energy_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc; Ng = grid.Ng

    # Steps 1-2: A-grid contravariant (numerical)
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

    # Step 3: Ghost extend (numerical)
    ua_ext, va_ext = extend_with_ghosts_vector(ua, va, grid)

    # Step 4: Build ArrayOp for face-position KE
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    o = Ng  # offset: real cell k maps to extended index k + Ng
    ud_c = const_wrap(unwrap(u_d)); vd_c = const_wrap(unwrap(v_d))
    ua_c = const_wrap(ua_ext); va_c = const_wrap(va_ext)

    # η-contribution: u_d · ṽ_η at UEdge faces
    # UEdge(i,j): va interpolated from cells (i-1,j) and (i,j)
    va_west = (wrap(va_c[p, i - 1 + o, j + o]) + wrap(va_c[p, i + o, j + o])) / 2
    # UEdge(i+1,j): va interpolated from cells (i,j) and (i+1,j)
    va_east = (wrap(va_c[p, i + o, j + o]) + wrap(va_c[p, i + 1 + o, j + o])) / 2
    ke_eta = (wrap(ud_c[p, i, j]) * va_west + wrap(ud_c[p, i + 1, j]) * va_east) / 2

    # ξ-contribution: v_d · ũ_ξ at VEdge faces
    # VEdge(i,j): ua interpolated from cells (i,j-1) and (i,j)
    ua_south = (wrap(ua_c[p, i + o, j - 1 + o]) + wrap(ua_c[p, i + o, j + o])) / 2
    # VEdge(i,j+1): ua interpolated from cells (i,j) and (i,j+1)
    ua_north = (wrap(ua_c[p, i + o, j + o]) + wrap(ua_c[p, i + o, j + 1 + o])) / 2
    ke_xi = (wrap(vd_c[p, i, j]) * ua_south + wrap(vd_c[p, i, j + 1]) * ua_north) / 2

    expr = (ke_xi + ke_eta) / 2
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end

# =============================================================================
# Cell-center KE (alternative formulation)
# =============================================================================

"""
    fv_kinetic_energy_cell!(ke, u_d, v_d, grid)

Compute cell-center kinetic energy from D-grid normalized covariant winds
using simple cell-center averaging.

D-grid winds are averaged to cell centers, then KE is computed using
the non-orthogonal formula with the center sin/cos(α) (super-grid position 5).

The averaging uses the 4 surrounding D-grid wind values:
    u_center = (u_d[i, j] + u_d[i+1, j]) / 2      (η-component, from UEdge)
    v_center = (v_d[i, j] + v_d[i, j+1]) / 2      (ξ-component, from VEdge)

Then: KE = (u² + v² - 2uv·cos α) / (2·sin²α)

where u and v are the normalized covariant components (V·ê_η and V·ê_ξ).

Note: This cell-center formulation can trigger the Hollingsworth-Kallberg
instability in dynamical core simulations. Prefer `fv_kinetic_energy!` for
production use.
"""
function fv_kinetic_energy_cell!(ke, u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Average D-grid winds to cell center
        # u_d (η-component) at west and east edges of cell
        u_avg = 0.5 * (u_d[p, i, j] + u_d[p, i + 1, j])
        # v_d (ξ-component) at south and north edges of cell
        v_avg = 0.5 * (v_d[p, i, j] + v_d[p, i, j + 1])

        # Get center sin/cos(α) from super-grid position 5
        sin_a = grid.sin_sg[p, i, j, 5]
        cos_a = grid.cos_sg[p, i, j, 5]
        sin2 = sin_a^2

        # KE = (u² + v² - 2uv·cos α) / (2·sin²α)
        ke[p, i, j] = (u_avg^2 + v_avg^2 - 2 * u_avg * v_avg * cos_a) / (2 * sin2)
    end

    return ke
end

"""
    fv_kinetic_energy_cell(u_d, v_d, grid)

Allocating version of `fv_kinetic_energy_cell!`. Returns ke array (6, Nc, Nc).

Note: This cell-center formulation can trigger the Hollingsworth-Kallberg
instability. Prefer `fv_kinetic_energy` for production use.
"""
function fv_kinetic_energy_cell(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    T = promote_type(eltype(u_d), eltype(v_d))
    ke = zeros(T, 6, Nc, Nc)
    fv_kinetic_energy_cell!(ke, u_d, v_d, grid)
    return ke
end

"""
    fv_kinetic_energy_cell_arrayop(u_d, v_d, grid)

ArrayOp version of the cell-center kinetic energy operator for the symbolic
discretization pipeline. Returns an ArrayOp{SymReal} of shape [6, Nc, Nc].

Note: This cell-center formulation can trigger the Hollingsworth-Kallberg
instability. Prefer `fv_kinetic_energy_arrayop` for production use.
"""
function fv_kinetic_energy_cell_arrayop(u_d, v_d, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    u_c = const_wrap(unwrap(u_d)); v_c = const_wrap(unwrap(v_d))
    sin_c = const_wrap(grid.sin_sg[:, :, :, 5])
    cos_c = const_wrap(grid.cos_sg[:, :, :, 5])

    # Average D-grid to cell center
    u_avg = (wrap(u_c[p, i, j]) + wrap(u_c[p, i + 1, j])) / 2
    v_avg = (wrap(v_c[p, i, j]) + wrap(v_c[p, i, j + 1])) / 2

    sin_a = wrap(sin_c[p, i, j])
    cos_a = wrap(cos_c[p, i, j])

    expr = (u_avg^2 + v_avg^2 - 2 * u_avg * v_avg * cos_a) / (2 * sin_a^2)
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end
