@testsnippet ShallowWaterSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Williamson test case: grid area" setup = [ShallowWaterSetup] tags = [:shallow_water] begin
    R = 6.371e6
    for Nc in [8, 16]
        grid = CubedSphereGrid(Nc; R = R)
        @test grid.Nc == Nc
        @test isapprox(total_area(grid), 4pi * R^2; rtol = 1.0e-10)
    end
end

@testitem "Williamson case 1: cosine bell initial condition" setup = [ShallowWaterSetup] tags = [:shallow_water] begin
    R = 6.371e6
    Nc = 16
    grid = CubedSphereGrid(Nc; R = R)

    # Cosine bell centered at (lon0, lat0) = (0, 0)
    r0 = R / 3  # bell radius
    h0 = 1000.0  # bell height
    lon0, lat0 = 0.0, 0.0

    h = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        lo = grid.lon[p, i, j]
        la = grid.lat[p, i, j]
        r = R * acos(clamp(sin(lat0) * sin(la) + cos(lat0) * cos(la) * cos(lo - lon0), -1.0, 1.0))
        h[p, i, j] = r < r0 ? (h0 / 2) * (1 + cos(pi * r / r0)) : 0.0
    end

    @test size(h) == (6, Nc, Nc)
    @test maximum(h) <= h0
    @test maximum(h) > h0 * 0.5  # Bell peak should be significant

    # Mass conservation check: total mass = ∫ h dA
    total_mass = sum(h[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    @test total_mass > 0
end
