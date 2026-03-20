"""
Finite-volume Laplacian with proper metric corrections including
the full off-diagonal g^{ξη} cross-derivative terms.

The full covariant Laplacian on the cubed sphere is:
    ∇²φ = (1/J)[∂/∂ξ(J g^{ξξ} ∂φ/∂ξ + J g^{ξη} ∂φ/∂η) +
                 ∂/∂η(J g^{ξη} ∂φ/∂ξ + J g^{ηη} ∂φ/∂η)]

The orthogonal part (g^{ξξ}, g^{ηη}) is computed using physical
center-to-center distances and edge lengths.

The cross-metric terms expand to:
    (1/J)[∂/∂ξ(Jg^{ξη})·∂φ/∂η + Jg^{ξη}·∂²φ/∂ξ∂η
        + ∂/∂η(Jg^{ξη})·∂φ/∂ξ + Jg^{ξη}·∂²φ/∂ξ∂η]

  = 2·g^{ξη}·∂²φ/(∂ξ∂η)
  + (1/J)·∂(Jg^{ξη})/∂ξ · ∂φ/∂η
  + (1/J)·∂(Jg^{ξη})/∂η · ∂φ/∂ξ
"""

"""
    fv_laplacian(phi, grid)

ArrayOp for the Laplacian at interior cells [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
Includes the complete off-diagonal metric correction for non-orthogonal grids,
including the gradient-of-metric terms ∂(Jg^{ξη})/∂ξ and ∂(Jg^{ξη})/∂η.
"""
function fv_laplacian(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi)); A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    dist_xi_c = const_wrap(grid.dist_xi); dist_eta_c = const_wrap(grid.dist_eta)
    J_c = const_wrap(grid.J); gxe_c = const_wrap(grid.ginv_ξη)
    dJgxe_dξ_c = const_wrap(grid.dJgxe_dξ)
    dJgxe_dη_c = const_wrap(grid.dJgxe_dη)

    dξ = grid.dξ; dη = grid.dη

    # Orthogonal part: 5-point stencil using physical distances and edge lengths
    orthogonal = wrap(1 / A_c[p, i + 1, j + 1]) * (
        (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_xi_c[p, i + 1, j + 1]) * wrap(dx_c[p, i + 2, j + 1]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i, j + 1]))     / wrap(dist_xi_c[p, i, j + 1])     * wrap(dx_c[p, i + 1, j + 1]) +
        (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_eta_c[p, i + 1, j + 1]) * wrap(dy_c[p, i + 1, j + 2]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i + 1, j]))     / wrap(dist_eta_c[p, i + 1, j])     * wrap(dy_c[p, i + 1, j + 1]))

    # Mixed derivative ∂²φ/(∂ξ∂η) via 4-point cross stencil
    d2phi_dxideta = (wrap(phi_c[p, i + 2, j + 2]) - wrap(phi_c[p, i + 2, j]) -
                     wrap(phi_c[p, i, j + 2]) + wrap(phi_c[p, i, j])) / (4 * dξ * dη)

    # First derivatives ∂φ/∂ξ and ∂φ/∂η via centered differences
    dphi_dxi  = (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i, j + 1])) / (2 * dξ)
    dphi_deta = (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j])) / (2 * dη)

    # Full cross-metric correction:
    # 2·g^{ξη}·∂²φ/(∂ξ∂η) + (1/J)·∂(Jg^{ξη})/∂ξ·∂φ/∂η + (1/J)·∂(Jg^{ξη})/∂η·∂φ/∂ξ
    cross_term = 2 * wrap(gxe_c[p, i + 1, j + 1]) * d2phi_dxideta +
        wrap(1 / J_c[p, i + 1, j + 1]) * (
            wrap(dJgxe_dξ_c[p, i + 1, j + 1]) * dphi_deta +
            wrap(dJgxe_dη_c[p, i + 1, j + 1]) * dphi_dxi)

    expr = orthogonal + cross_term
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
