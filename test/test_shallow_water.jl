@testsnippet ShallowWaterSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Williamson test case placeholder" setup=[ShallowWaterSetup] tags=[:shallow_water] begin
    R = 6.371e6
    for Nc in [8, 16]
        grid = CubedSphereGrid(Nc; R=R)
        @test grid.Nc == Nc
        @test isapprox(total_area(grid), 4pi * R^2; rtol=1e-10)
    end
end
