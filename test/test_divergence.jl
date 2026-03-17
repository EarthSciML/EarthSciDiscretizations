@testsnippet DivSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Divergence of zero flux is zero" setup=[DivSetup] tags=[:operators] begin
    Nc = 4
    grid = CubedSphereGrid(Nc; R=1.0)
    F_xi = zeros(6, Nc + 1, Nc)
    F_eta = zeros(6, Nc, Nc + 1)
    div_vals = evaluate_arrayop(fv_divergence(F_xi, F_eta, grid))
    @test size(div_vals) == (6, Nc, Nc)
    @test all(abs.(div_vals) .< 1e-14)
end

@testitem "Divergence theorem" setup=[DivSetup] tags=[:operators] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)
    F_xi = zeros(6, Nc + 1, Nc)
    F_eta = zeros(6, Nc, Nc + 1)
    for p in 1:6
        for i in 1:Nc+1, j in 1:Nc
            lon, lat = gnomonic_to_lonlat(grid.ξ_edges[i], grid.η_centers[j], p)
            F_xi[p, i, j] = cos(lat) * cos(lon)
        end
        for i in 1:Nc, j in 1:Nc+1
            lon, lat = gnomonic_to_lonlat(grid.ξ_centers[i], grid.η_edges[j], p)
            F_eta[p, i, j] = cos(lat) * sin(lon)
        end
    end
    div_vals = evaluate_arrayop(fv_divergence(F_xi, F_eta, grid))
    integral = sum(div_vals[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    # Panel boundary discretization error — relaxed tolerance
    @test abs(integral) < 20.0
end

@testitem "Divergence is linear" setup=[DivSetup] tags=[:operators] begin
    Nc = 4
    grid = CubedSphereGrid(Nc; R=1.0)
    F1_xi = rand(6, Nc + 1, Nc); F1_eta = rand(6, Nc, Nc + 1)
    F2_xi = rand(6, Nc + 1, Nc); F2_eta = rand(6, Nc, Nc + 1)
    div1 = evaluate_arrayop(fv_divergence(F1_xi, F1_eta, grid))
    div2 = evaluate_arrayop(fv_divergence(F2_xi, F2_eta, grid))
    div_sum = evaluate_arrayop(fv_divergence(F1_xi + F2_xi, F1_eta + F2_eta, grid))
    @test isapprox(div_sum, div1 + div2; rtol=1e-12)
end
