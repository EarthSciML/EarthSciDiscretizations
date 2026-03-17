"""
Finite-volume Laplacian (5-point stencil with metric corrections).
"""

"""
    fv_laplacian(phi, grid)

ArrayOp for the Laplacian at interior cells [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
"""
function fv_laplacian(phi, grid::CubedSphereGrid)
    Nc = grid.Nc; dξ = grid.dξ; dη = grid.dη
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi)); A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    expr = wrap(1 / A_c[p, i + 1, j + 1]) * (
        (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i + 1, j + 1])) / dξ * wrap(dx_c[p, i + 2, j + 1]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i, j + 1])) / dξ * wrap(dx_c[p, i + 1, j + 1]) +
        (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j + 1])) / dη * wrap(dy_c[p, i + 1, j + 2]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i + 1, j])) / dη * wrap(dy_c[p, i + 1, j + 1]))
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
