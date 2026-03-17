@testsnippet LapSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Laplacian of constant is zero" setup=[LapSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)
    phi = fill(7.0, 6, Nc, Nc)
    lap = evaluate_arrayop(fv_laplacian(phi, grid))
    @test size(lap) == (6, Nc - 2, Nc - 2)
    @test all(abs.(lap) .< 1e-10)
end

@testitem "Laplacian convergence" setup=[LapSetup] tags=[:operators] begin
    errors = Float64[]
    for Nc in [8, 16, 32]
        grid = CubedSphereGrid(Nc; R=1.0)
        phi = zeros(6, Nc, Nc)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
        end
        lap = evaluate_arrayop(fv_laplacian(phi, grid))
        push!(errors, maximum(abs.(lap)))
    end
    # Error should decrease (or at least not grow) with resolution
    for k in 2:length(errors)
        @test errors[k] < errors[k-1] * 1.5
    end
end

@testitem "Laplacian is linear" setup=[LapSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)
    a = rand(6, Nc, Nc); b = rand(6, Nc, Nc); α = 3.5
    la = evaluate_arrayop(fv_laplacian(a, grid))
    lb = evaluate_arrayop(fv_laplacian(b, grid))
    lab = evaluate_arrayop(fv_laplacian(a + α * b, grid))
    @test isapprox(lab, la + α * lb; rtol=1e-10)
end
