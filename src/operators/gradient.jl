"""
Finite-volume gradient operator on the cubed-sphere C-grid.
"""

"""
    fv_gradient_xi(phi, grid)

ArrayOp for the ξ-gradient at interior edges [6, Nc-1, Nc].
"""
function fv_gradient_xi(phi, grid::CubedSphereGrid)
    Nc = grid.Nc; dξ = grid.dξ
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi))
    expr = (wrap(phi_c[p, i + 1, j]) - wrap(phi_c[p, i, j])) / dξ
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 1), j => 1:1:Nc))
end

"""
    fv_gradient_eta(phi, grid)

ArrayOp for the η-gradient at interior edges [6, Nc, Nc-1].
"""
function fv_gradient_eta(phi, grid::CubedSphereGrid)
    Nc = grid.Nc; dη = grid.dη
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi))
    expr = (wrap(phi_c[p, i, j + 1]) - wrap(phi_c[p, i, j])) / dη
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc - 1)))
end
