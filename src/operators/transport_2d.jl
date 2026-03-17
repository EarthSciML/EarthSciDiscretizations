"""
2D transport operator (Lax-Friedrichs upwind in both directions).
"""

"""
    transport_2d(q, courant_xi, courant_eta, grid)

ArrayOp for 2D transport tendency [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
"""
function transport_2d(q, courant_xi, courant_eta, grid::CubedSphereGrid)
    Nc = grid.Nc; dξ = grid.dξ; dη = grid.dη
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    q_c = const_wrap(unwrap(q)); cxi = const_wrap(unwrap(courant_xi)); ceta = const_wrap(unwrap(courant_eta))

    F_e_a = wrap(cxi[p, i + 2, j + 1]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 2, j + 1])) / 2
    F_e_d = abs(wrap(cxi[p, i + 2, j + 1])) * (wrap(q_c[p, i + 2, j + 1]) - wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_a = wrap(cxi[p, i + 1, j + 1]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i + 1, j + 1])) / 2
    F_w_d = abs(wrap(cxi[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i, j + 1])) / 2
    G_n_a = wrap(ceta[p, i + 1, j + 2]) * (wrap(q_c[p, i + 1, j + 1]) + wrap(q_c[p, i + 1, j + 2])) / 2
    G_n_d = abs(wrap(ceta[p, i + 1, j + 2])) * (wrap(q_c[p, i + 1, j + 2]) - wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_a = wrap(ceta[p, i + 1, j + 1]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 1, j + 1])) / 2
    G_s_d = abs(wrap(ceta[p, i + 1, j + 1])) * (wrap(q_c[p, i + 1, j + 1]) - wrap(q_c[p, i + 1, j])) / 2

    expr = -((F_e_a - F_e_d) - (F_w_a - F_w_d)) / dξ - ((G_n_a - G_n_d) - (G_s_a - G_s_d)) / dη
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
