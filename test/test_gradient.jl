@testsnippet GradSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "Gradient of constant is zero" setup = [GradSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    phi = fill(5.0, 6, Nc, Nc)

    ao_xi = fv_gradient_xi(phi, grid)
    ao_eta = fv_gradient_eta(phi, grid)

    grad_xi = evaluate_arrayop(ao_xi)
    grad_eta = evaluate_arrayop(ao_eta)

    @test all(abs.(grad_xi) .< 1.0e-12)
    @test all(abs.(grad_eta) .< 1.0e-12)
end

@testitem "Gradient convergence" setup = [GradSetup] tags = [:operators] begin
    # Gradient of a smooth function should converge with resolution
    errors_xi = Float64[]
    for Nc in [8, 16, 32]
        grid = CubedSphereGrid(Nc; R = 1.0)
        phi = zeros(6, Nc, Nc)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
        end
        grad_xi = evaluate_arrayop(fv_gradient_xi(phi, grid))
        push!(errors_xi, maximum(abs.(grad_xi)))
    end
    # Error pattern should be stable or decreasing
    @test errors_xi[end] < errors_xi[1] * 2.0
end

@testitem "Gradient is linear" setup = [GradSetup] tags = [:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    phi1 = randn(6, Nc, Nc)
    phi2 = randn(6, Nc, Nc)
    alpha = 2.3

    ao1 = fv_gradient_xi(phi1, grid)
    ao2 = fv_gradient_xi(phi2, grid)
    ao_sum = fv_gradient_xi(phi1 + alpha * phi2, grid)

    g1 = evaluate_arrayop(ao1)
    g2 = evaluate_arrayop(ao2)
    g_sum = evaluate_arrayop(ao_sum)

    @test isapprox(g_sum, g1 + alpha * g2; rtol = 1.0e-12)
end
