"""
1D flux-form transport (Lax-Friedrichs upwind).
"""

"""
    flux_1d(q, courant, grid, dim)

ArrayOp for 1D transport tendency. Index (p,i,j) maps to physical cell (p, i+1, j)
for `:xi` or (p, i, j+1) for `:eta`.
"""
function flux_1d(q, courant, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc; idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    q_c = const_wrap(unwrap(q)); c_c = const_wrap(unwrap(courant))
    if dim == :xi
        dξ = grid.dξ
        F_e_a = wrap(c_c[p, i + 2, j]) * (wrap(q_c[p, i + 1, j]) + wrap(q_c[p, i + 2, j])) / 2
        F_e_d = abs(wrap(c_c[p, i + 2, j])) * (wrap(q_c[p, i + 2, j]) - wrap(q_c[p, i + 1, j])) / 2
        F_w_a = wrap(c_c[p, i + 1, j]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i + 1, j])) / 2
        F_w_d = abs(wrap(c_c[p, i + 1, j])) * (wrap(q_c[p, i + 1, j]) - wrap(q_c[p, i, j])) / 2
        expr = -((F_e_a - F_e_d) - (F_w_a - F_w_d)) / dξ
        ranges = Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:Nc)
    else
        dη = grid.dη
        F_n_a = wrap(c_c[p, i, j + 2]) * (wrap(q_c[p, i, j + 1]) + wrap(q_c[p, i, j + 2])) / 2
        F_n_d = abs(wrap(c_c[p, i, j + 2])) * (wrap(q_c[p, i, j + 2]) - wrap(q_c[p, i, j + 1])) / 2
        F_s_a = wrap(c_c[p, i, j + 1]) * (wrap(q_c[p, i, j]) + wrap(q_c[p, i, j + 1])) / 2
        F_s_d = abs(wrap(c_c[p, i, j + 1])) * (wrap(q_c[p, i, j + 1]) - wrap(q_c[p, i, j])) / 2
        expr = -((F_n_a - F_n_d) - (F_s_a - F_s_d)) / dη
        ranges = Dict(p => 1:1:6, i => 1:1:Nc, j => 1:1:(Nc - 2))
    end
    return make_arrayop(idx, unwrap(expr), ranges)
end
