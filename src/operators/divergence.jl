"""
Finite-volume divergence operator on the cubed-sphere C-grid.

    div(F) = (1/A) * [F_ξ·dx|_{i+1} - F_ξ·dx|_{i} + F_η·dy|_{j+1} - F_η·dy|_{j}]
"""

"""
    fv_divergence(F_xi, F_eta, grid)

ArrayOp for the finite-volume divergence at CellCenter [6, Nc, Nc].
"""
function fv_divergence(F_xi, F_eta, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    A_c = const_wrap(grid.area); dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    Fxi_c = const_wrap(unwrap(F_xi)); Feta_c = const_wrap(unwrap(F_eta))
    expr = wrap(1 / A_c[p, i, j]) * (
        wrap(Fxi_c[p, i + 1, j]) * wrap(dx_c[p, i + 1, j]) - wrap(Fxi_c[p, i, j]) * wrap(dx_c[p, i, j]) +
            wrap(Feta_c[p, i, j + 1]) * wrap(dy_c[p, i, j + 1]) - wrap(Feta_c[p, i, j]) * wrap(dy_c[p, i, j])
    )
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:Nc))
end
