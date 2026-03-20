"""
Full covariant Laplacian on the cubed sphere:
    ∇²φ = (1/J)[∂/∂ξ(J g^{ξξ} ∂φ/∂ξ + J g^{ξη} ∂φ/∂η) +
                 ∂/∂η(J g^{ξη} ∂φ/∂ξ + J g^{ηη} ∂φ/∂η)]

Expanding all terms:
  Orthogonal:
    g^{ξξ}·∂²φ/∂ξ² + (1/J)·∂(Jg^{ξξ})/∂ξ·∂φ/∂ξ
  + g^{ηη}·∂²φ/∂η² + (1/J)·∂(Jg^{ηη})/∂η·∂φ/∂η
  Cross-metric:
    2·g^{ξη}·∂²φ/(∂ξ∂η)
  + (1/J)·∂(Jg^{ξη})/∂ξ · ∂φ/∂η
  + (1/J)·∂(Jg^{ξη})/∂η · ∂φ/∂ξ
"""

"""
    fv_laplacian(phi, grid)

ArrayOp for the Laplacian at interior cells [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
Uses the full covariant form with explicit inverse metric tensor
components g^{ξξ}, g^{ηη}, g^{ξη and their derivatives, ensuring
correctness on the non-orthogonal cubed-sphere grid.
"""
function fv_laplacian(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi))
    J_c = const_wrap(grid.J)
    ginv_ξξ_c = const_wrap(grid.ginv_ξξ)
    ginv_ηη_c = const_wrap(grid.ginv_ηη)
    gxe_c = const_wrap(grid.ginv_ξη)
    dJgxx_dξ_c = const_wrap(grid.dJgxx_dξ)
    dJgyy_dη_c = const_wrap(grid.dJgyy_dη)
    dJgxe_dξ_c = const_wrap(grid.dJgxe_dξ)
    dJgxe_dη_c = const_wrap(grid.dJgxe_dη)

    dξ = grid.dξ; dη = grid.dη

    # Second derivatives in computational coordinates (centered)
    d2phi_dξ2 = (wrap(phi_c[p, i + 2, j + 1]) - 2 * wrap(phi_c[p, i + 1, j + 1]) +
                 wrap(phi_c[p, i, j + 1])) / dξ^2
    d2phi_dη2 = (wrap(phi_c[p, i + 1, j + 2]) - 2 * wrap(phi_c[p, i + 1, j + 1]) +
                 wrap(phi_c[p, i + 1, j])) / dη^2

    # First derivatives (centered)
    dphi_dξ = (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i, j + 1])) / (2 * dξ)
    dphi_dη = (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j])) / (2 * dη)

    # Orthogonal part:
    # (1/J)∂/∂ξ(Jg^{ξξ}∂φ/∂ξ) = g^{ξξ}·∂²φ/∂ξ² + (1/J)·∂(Jg^{ξξ})/∂ξ·∂φ/∂ξ
    # (1/J)∂/∂η(Jg^{ηη}∂φ/∂η) = g^{ηη}·∂²φ/∂η² + (1/J)·∂(Jg^{ηη})/∂η·∂φ/∂η
    orthogonal = wrap(ginv_ξξ_c[p, i + 1, j + 1]) * d2phi_dξ2 +
                 wrap(ginv_ηη_c[p, i + 1, j + 1]) * d2phi_dη2 +
        wrap(1 / J_c[p, i + 1, j + 1]) * (
            wrap(dJgxx_dξ_c[p, i + 1, j + 1]) * dphi_dξ +
            wrap(dJgyy_dη_c[p, i + 1, j + 1]) * dphi_dη)

    # Mixed derivative ∂²φ/(∂ξ∂η) via 4-point cross stencil
    d2phi_dξdη = (wrap(phi_c[p, i + 2, j + 2]) - wrap(phi_c[p, i + 2, j]) -
                  wrap(phi_c[p, i, j + 2]) + wrap(phi_c[p, i, j])) / (4 * dξ * dη)

    # Cross-metric correction:
    # 2·g^{ξη}·∂²φ/(∂ξ∂η) + (1/J)·∂(Jg^{ξη})/∂ξ·∂φ/∂η + (1/J)·∂(Jg^{ξη})/∂η·∂φ/∂ξ
    cross_term = 2 * wrap(gxe_c[p, i + 1, j + 1]) * d2phi_dξdη +
        wrap(1 / J_c[p, i + 1, j + 1]) * (
            wrap(dJgxe_dξ_c[p, i + 1, j + 1]) * dphi_dη +
            wrap(dJgxe_dη_c[p, i + 1, j + 1]) * dphi_dξ)

    expr = orthogonal + cross_term
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
