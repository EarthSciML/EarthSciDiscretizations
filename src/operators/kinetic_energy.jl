"""
FV3 kinetic energy computation on the non-orthogonal cubed-sphere grid.

On a non-orthogonal grid with angle α between unit coordinate vectors ê_ξ and ê_η,
the kinetic energy is computed from normalized covariant components (u = V·ê_ξ,
v = V·ê_η):

    KE = (u² + v² - 2uv·cos α) / (2·sin²α)

This is equivalent to the mixed form KE = (u·ũ + v·ṽ) / 2, where ũ, ṽ are
the contravariant components.

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Eq. 3.5.

## FV3 upstream-biased KE (Section 5)

FV3 computes KE in an upstream-biased manner to avoid the Hollingsworth-Kallberg
instability that arises from computing KE on the A-grid from D-grid winds:

    KE* = (1/2) · [X(u, ũ*_b) + Y(v, ṽ*_b)]

where X and Y are 1D transport operators and ũ*_b, ṽ*_b are advective winds
averaged to grid corners. This approach naturally damps the computational mode.

For now we implement the simpler cell-center KE which is sufficient for many
applications.
"""

"""
    fv_kinetic_energy!(ke, u_d, v_d, grid)

Compute cell-center kinetic energy from D-grid normalized covariant winds.

D-grid winds are averaged to cell centers, then KE is computed using
the non-orthogonal formula with the center sin/cos(α) (super-grid position 5).

The averaging uses the 4 surrounding D-grid wind values:
    u_center = (u_d[i, j] + u_d[i+1, j]) / 2      (η-component, from UEdge)
    v_center = (v_d[i, j] + v_d[i, j+1]) / 2      (ξ-component, from VEdge)

Then: KE = (u² + v² - 2uv·cos α) / (2·sin²α)

where u and v are the normalized covariant components (V·ê_η and V·ê_ξ).
"""
function fv_kinetic_energy!(ke, u_d, v_d, grid::CubedSphereGrid)
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
    fv_kinetic_energy(u_d, v_d, grid)

Allocating version of `fv_kinetic_energy!`. Returns ke array (6, Nc, Nc).
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

ArrayOp version of the FV3 kinetic energy operator for the symbolic
discretization pipeline. Returns an ArrayOp{SymReal} of shape [6, Nc, Nc].
"""
function fv_kinetic_energy_arrayop(u_d, v_d, grid::CubedSphereGrid)
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
