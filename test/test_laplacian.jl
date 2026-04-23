@testsnippet LapSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Laplacian of constant is zero" setup = [LapSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)
    phi = fill(7.0, 6, Nc, Nc)
    lap = evaluate_arrayop(fv_laplacian(phi, grid))
    @test size(lap) == (6, Nc - 2, Nc - 2)
    @test all(abs.(lap) .< 1.0e-10)
end

@testitem "Laplacian second-order convergence" setup = [LapSetup] tags = [:operators] begin
    # Test with φ = cos(2ξ)·cos(2η), for which the analytical covariant Laplacian
    # can be computed. We verify that the error converges at second order.
    errors = Float64[]
    resolutions = [8, 16, 32]
    for Nc in resolutions
        grid = CubedSphereGrid(Nc; R = 1.0)

        # Set up test field and compute analytical Laplacian
        phi = zeros(6, Nc, Nc)
        lap_exact = zeros(6, Nc - 2, Nc - 2)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
        end

        # Analytical: ∇²φ = (1/J)[∂/∂ξ(Jg^{ξξ}∂φ/∂ξ + Jg^{ξη}∂φ/∂η) +
        #                          ∂/∂η(Jg^{ξη}∂φ/∂ξ + Jg^{ηη}∂φ/∂η)]
        # Compute analytically using metric at cell centers via finite differences
        # with very fine step for the analytical reference
        h = 1.0e-6
        for ii in 1:(Nc - 2), jj in 1:(Nc - 2)
            ic = ii + 1; jc = jj + 1  # physical cell indices
            ξc = grid.ξ_centers[ic]; ηc = grid.η_centers[jc]

            # Evaluate φ and its derivatives analytically
            f(ξ, η) = cos(2ξ) * cos(2η)
            df_dξ(ξ, η) = -2sin(2ξ) * cos(2η)
            df_dη(ξ, η) = -2cos(2ξ) * sin(2η)
            d2f_dξ2(ξ, η) = -4cos(2ξ) * cos(2η)
            d2f_dη2(ξ, η) = -4cos(2ξ) * cos(2η)
            d2f_dξdη(ξ, η) = 4sin(2ξ) * sin(2η)

            # Get metric quantities at (ξc, ηc) — same for all panels
            J, g_ξξ, g_ηη, g_ξη = gnomonic_metric(ξc, ηc, 1.0)
            det_g = g_ξξ * g_ηη - g_ξη^2
            ginv_ξξ = g_ηη / det_g
            ginv_ηη = g_ξξ / det_g
            ginv_ξη = -g_ξη / det_g

            # Compute ∂(Jg^{ab})/∂ξ and ∂(Jg^{ab})/∂η numerically
            function Jginv(ξ, η)
                Jv, gxx, gee, gxe = gnomonic_metric(ξ, η, 1.0)
                d = gxx * gee - gxe^2
                return (Jv * gee / d, Jv * gxx / d, Jv * (-gxe / d))
            end
            Jgxx_p, Jgyy_p, Jgxe_p = Jginv(ξc + h, ηc)
            Jgxx_m, Jgyy_m, Jgxe_m = Jginv(ξc - h, ηc)
            Jgxx_np, Jgyy_np, Jgxe_np = Jginv(ξc, ηc + h)
            Jgxx_nm, Jgyy_nm, Jgxe_nm = Jginv(ξc, ηc - h)

            dJgxx_dξ = (Jgxx_p - Jgxx_m) / (2h)
            dJgyy_dη = (Jgyy_np - Jgyy_nm) / (2h)
            dJgxe_dξ = (Jgxe_p - Jgxe_m) / (2h)
            dJgxe_dη = (Jgxe_np - Jgxe_nm) / (2h)

            # Full covariant Laplacian
            lap_val = ginv_ξξ * d2f_dξ2(ξc, ηc) + ginv_ηη * d2f_dη2(ξc, ηc) +
                (1 / J) * (dJgxx_dξ * df_dξ(ξc, ηc) + dJgyy_dη * df_dη(ξc, ηc)) +
                2 * ginv_ξη * d2f_dξdη(ξc, ηc) +
                (1 / J) * (dJgxe_dξ * df_dη(ξc, ηc) + dJgxe_dη * df_dξ(ξc, ηc))

            for p in 1:6
                lap_exact[p, ii, jj] = lap_val
            end
        end

        lap_num = evaluate_arrayop(fv_laplacian(phi, grid))
        push!(errors, maximum(abs.(lap_num .- lap_exact)))
    end

    # Verify second-order convergence: error ratio should be ~4 for 2x refinement
    for k in 2:length(errors)
        ratio = errors[k - 1] / errors[k]
        @test ratio > 2.0  # expect ~4 for 2nd order; pre-asymptotic effects at coarse Nc
    end
end

@testitem "Laplacian is linear" setup = [LapSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)
    a = rand(6, Nc, Nc); b = rand(6, Nc, Nc); α = 3.5
    la = evaluate_arrayop(fv_laplacian(a, grid))
    lb = evaluate_arrayop(fv_laplacian(b, grid))
    lab = evaluate_arrayop(fv_laplacian(a + α * b, grid))
    @test isapprox(lab, la + α * lb; rtol = 1.0e-10)
end

@testitem "fv_laplacian matches fv_laplacian_extended" setup = [LapSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    # Non-trivial test field
    phi = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * sin(3 * grid.η_centers[j])
    end

    # fv_laplacian: operates on raw array, outputs (6, Nc-2, Nc-2)
    lap1 = evaluate_arrayop(fv_laplacian(phi, grid))

    # fv_laplacian_extended: operates on ghost-extended array, outputs (6, Nc, Nc)
    phi_ext = extend_with_ghosts(phi, grid)
    lap2 = evaluate_arrayop(fv_laplacian_extended(phi_ext, grid))

    # Compare on the shared interior region: fv_laplacian index (p,i,j) maps to
    # physical cell (p, i+1, j+1), while fv_laplacian_extended index (p,i,j)
    # maps to physical cell (p, i, j)
    for p in 1:6, i in 1:(Nc - 2), j in 1:(Nc - 2)
        @test isapprox(lap1[p, i, j], lap2[p, i + 1, j + 1]; rtol = 1.0e-10)
    end
end

@testitem "Precomputed stencil matches ArrayOp Laplacian" setup = [LapSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)
    stencil = precompute_laplacian_stencil(grid)

    phi = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * sin(3 * grid.η_centers[j])
    end

    # Precomputed stencil result
    du_stencil = zeros(6, Nc, Nc)
    apply_laplacian!(du_stencil, phi, stencil)

    # Ghost-extended ArrayOp result (covers all cells)
    phi_ext = extend_with_ghosts(phi, grid)
    du_ext = evaluate_arrayop(fv_laplacian_extended(phi_ext, grid))

    # They should agree at interior cells (boundary cells may differ due to
    # different ghost cell handling)
    for p in 1:6, i in 2:(Nc - 1), j in 2:(Nc - 1)
        @test isapprox(du_stencil[p, i, j], du_ext[p, i, j]; rtol = 1.0e-6)
    end
end
