"""
Finite-volume Laplacian with proper metric corrections including
the off-diagonal g^{ξη} cross-derivative term.

The full covariant Laplacian on the cubed sphere is:
    ∇²φ = (1/J)[∂/∂ξ(J g^{ξξ} ∂φ/∂ξ + J g^{ξη} ∂φ/∂η) +
                 ∂/∂η(J g^{ξη} ∂φ/∂ξ + J g^{ηη} ∂φ/∂η)]

The orthogonal part (g^{ξξ}, g^{ηη}) is computed using physical
center-to-center distances and edge lengths. The cross-metric
g^{ξη} correction includes:
  1. The mixed second derivative: 2 g^{ξη} ∂²φ/(∂ξ∂η)
  2. First-derivative corrections from spatial variation of J·g^{ξη}:
     (1/J)[∂(J g^{ξη})/∂ξ · ∂φ/∂η + ∂(J g^{ξη})/∂η · ∂φ/∂ξ]

    ∇²φ ≈ (1/A) Σ_edges [(φ_neighbor - φ_center) / dist_phys * edge_length]
         + 2 g^{ξη} ∂²φ/(∂ξ∂η)
         + (1/J)[∂(J g^{ξη})/∂ξ · ∂φ/∂η + ∂(J g^{ξη})/∂η · ∂φ/∂ξ]
"""

"""
    fv_laplacian(phi, grid)

ArrayOp for the Laplacian at interior cells [6, Nc-2, Nc-2].
Index (p,i,j) maps to physical cell (p, i+1, j+1).
Includes the off-diagonal metric correction for non-orthogonal grids.
"""
function fv_laplacian(phi, grid::CubedSphereGrid)
    Nc = grid.Nc
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    phi_c = const_wrap(unwrap(phi)); A_c = const_wrap(grid.area)
    dx_c = const_wrap(grid.dx); dy_c = const_wrap(grid.dy)
    dist_xi_c = const_wrap(grid.dist_xi); dist_eta_c = const_wrap(grid.dist_eta)
    J_c = const_wrap(grid.J); gxe_c = const_wrap(grid.ginv_ξη)

    dξ = grid.dξ; dη = grid.dη

    # Precompute J * g^{ξη} product for first-derivative corrections
    Jgxe = grid.J .* grid.ginv_ξη
    Jgxe_c = const_wrap(Jgxe)

    # Orthogonal part: 5-point stencil using physical distances and edge lengths
    # East face: gradient uses dist_xi between cells (i+1,j+1) and (i+2,j+1), edge length dx at ξ-edge i+2
    # West face: gradient uses dist_xi between cells (i,j+1) and (i+1,j+1), edge length dx at ξ-edge i+1
    # North face: gradient uses dist_eta between cells (i+1,j+1) and (i+1,j+2), edge length dy at η-edge j+2
    # South face: gradient uses dist_eta between cells (i+1,j) and (i+1,j+1), edge length dy at η-edge j+1
    orthogonal = wrap(1 / A_c[p, i + 1, j + 1]) * (
        (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_xi_c[p, i + 1, j + 1]) * wrap(dx_c[p, i + 2, j + 1]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i, j + 1]))     / wrap(dist_xi_c[p, i, j + 1])     * wrap(dx_c[p, i + 1, j + 1]) +
        (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j + 1])) / wrap(dist_eta_c[p, i + 1, j + 1]) * wrap(dy_c[p, i + 1, j + 2]) -
        (wrap(phi_c[p, i + 1, j + 1]) - wrap(phi_c[p, i + 1, j]))     / wrap(dist_eta_c[p, i + 1, j])     * wrap(dy_c[p, i + 1, j + 1]))

    # Cross-metric correction: 2 * g^{ξη} * ∂²φ/(∂ξ∂η)
    # From the covariant Laplacian: (1/J)[∂/∂ξ(J g^{ξη} ∂φ/∂η) + ∂/∂η(J g^{ξη} ∂φ/∂ξ)]
    # Under constant J·g^{ξη} assumption, both terms reduce to J·g^{ξη}·∂²φ/(∂ξ∂η),
    # and (1/J) times their sum gives 2·g^{ξη}·∂²φ/(∂ξ∂η).
    cross_term = 2 * wrap(gxe_c[p, i + 1, j + 1]) *
        (wrap(phi_c[p, i + 2, j + 2]) - wrap(phi_c[p, i + 2, j]) -
         wrap(phi_c[p, i, j + 2]) + wrap(phi_c[p, i, j])) / (4 * dξ * dη)

    # First-derivative corrections from spatial variation of J·g^{ξη}:
    # (1/J)[∂(J g^{ξη})/∂ξ · ∂φ/∂η + ∂(J g^{ξη})/∂η · ∂φ/∂ξ]
    dJgxe_dξ = (wrap(Jgxe_c[p, i + 2, j + 1]) - wrap(Jgxe_c[p, i, j + 1])) / (2 * dξ)
    dJgxe_dη = (wrap(Jgxe_c[p, i + 1, j + 2]) - wrap(Jgxe_c[p, i + 1, j])) / (2 * dη)
    dphi_dξ = (wrap(phi_c[p, i + 2, j + 1]) - wrap(phi_c[p, i, j + 1])) / (2 * dξ)
    dphi_dη = (wrap(phi_c[p, i + 1, j + 2]) - wrap(phi_c[p, i + 1, j])) / (2 * dη)
    first_deriv_correction = wrap(1 / J_c[p, i + 1, j + 1]) * (dJgxe_dξ * dphi_dη + dJgxe_dη * dphi_dξ)

    expr = orthogonal + cross_term + first_deriv_correction
    return make_arrayop(idx, unwrap(expr), Dict(p => 1:1:6, i => 1:1:(Nc - 2), j => 1:1:(Nc - 2)))
end
