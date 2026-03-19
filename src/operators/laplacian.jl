"""
Finite-volume Laplacian with proper metric corrections.

Uses physical center-to-center distances for the normal gradient
at each cell face, multiplied by the physical edge length, and
normalized by cell area. This correctly handles the non-uniform
metric of the gnomonic cubed-sphere projection.

    ∇²φ ≈ (1/A) Σ_edges [(φ_neighbor - φ_center) / dist_phys * edge_length]
"""

"""
    fv_laplacian(phi, grid)

ArrayOp for the Laplacian at interior cells [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
"""
function fv_laplacian(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi)); A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    dist_xi_c = const_wrap(grid.dist_xi); dist_eta_c = const_wrap(grid.dist_eta)

    # East face: gradient uses dist_xi between cells (i+1,j+1) and (i+2,j+1), edge length dx at ξ-edge i+2
    # West face: gradient uses dist_xi between cells (i,j+1) and (i+1,j+1), edge length dx at ξ-edge i+1
    # North face: gradient uses dist_eta between cells (i+1,j+1) and (i+1,j+2), edge length dy at η-edge j+2
    # South face: gradient uses dist_eta between cells (i+1,j) and (i+1,j+1), edge length dy at η-edge j+1
    expr = wrap(1 / A_c[p, i + 1, j + 1]) * (
        (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_xi_c[p, i + 1, j + 1]) * wrap(dx_c[p, i + 2, j + 1]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i, j + 1]))     / wrap(dist_xi_c[p, i, j + 1])     * wrap(dx_c[p, i + 1, j + 1]) +
        (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_eta_c[p, i + 1, j + 1]) * wrap(dy_c[p, i + 1, j + 2]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i + 1, j]))     / wrap(dist_eta_c[p, i + 1, j])     * wrap(dy_c[p, i + 1, j + 1]))
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
