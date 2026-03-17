@testsnippet GradSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "Gradient of constant is zero" setup=[GradSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    phi = fill(5.0, 6, Nc, Nc)

    ao_xi = fv_gradient_xi(phi, grid)
    ao_eta = fv_gradient_eta(phi, grid)

    grad_xi = evaluate_arrayop(ao_xi)
    grad_eta = evaluate_arrayop(ao_eta)

    @test all(abs.(grad_xi) .< 1e-12)
    @test all(abs.(grad_eta) .< 1e-12)
end

@testitem "Gradient of linear xi-field is exact" setup=[GradSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    slope = 2.5

    # phi = slope * xi_center for each cell
    phi = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = slope * grid.ξ_centers[i]
    end

    ao_xi = fv_gradient_xi(phi, grid)
    grad_xi = evaluate_arrayop(ao_xi)

    # For a linear function in xi, the xi-gradient should be exactly the slope
    @test isapprox(grad_xi, fill(slope, size(grad_xi)...); rtol=1e-12)
end

@testitem "Gradient of linear eta-field is exact" setup=[GradSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    slope = -3.7

    # phi = slope * eta_center for each cell
    phi = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = slope * grid.η_centers[j]
    end

    ao_eta = fv_gradient_eta(phi, grid)
    grad_eta = evaluate_arrayop(ao_eta)

    # For a linear function in eta, the eta-gradient should be exactly the slope
    @test isapprox(grad_eta, fill(slope, size(grad_eta)...); rtol=1e-12)
end

@testitem "Gradient is linear" setup=[GradSetup] tags=[:operators] begin
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

    @test isapprox(g_sum, g1 + alpha * g2; rtol=1e-12)
end
