"""
Finite-volume gradient operator on the cubed-sphere grid.

Gradients are computed using physical center-to-center distances,
providing proper metric correction on the cubed sphere.
"""

"""
    fv_gradient_xi(phi, grid)

ArrayOp for the ξ-gradient at interior edges [6, Nc-1, Nc].
Uses physical center-to-center distance for metric correctness.
"""
function fv_gradient_xi(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi))
    dist_c = const_wrap(grid.dist_xi)
    expr = (wrap(phi_c[p, i + 1, j]) - wrap(phi_c[p, i, j])) / wrap(dist_c[p, i, j])
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 1), j => 1:1:Nc))
end

"""
    fv_gradient_eta(phi, grid)

ArrayOp for the η-gradient at interior edges [6, Nc, Nc-1].
Uses physical center-to-center distance for metric correctness.
"""
function fv_gradient_eta(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi))
    dist_c = const_wrap(grid.dist_eta)
    expr = (wrap(phi_c[p, i, j + 1]) - wrap(phi_c[p, i, j])) / wrap(dist_c[p, i, j])
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc - 1)))
end
