@testsnippet FV3Setup begin
    using Test
    using EarthSciDiscretizations
    using LinearAlgebra: norm, dot
end

# =============================================================================
# Super-grid (sin_sg, cos_sg) tests
# =============================================================================

@testitem "sin_sg: identity sin²α + cos²α = 1" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(16; R = 1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            s = grid.sin_sg[p, i, j, pos]
            c = grid.cos_sg[p, i, j, pos]
            @test isapprox(s^2 + c^2, 1.0; atol = 1.0e-14)
        end
    end
end

@testitem "sin_sg: positive sin(α) everywhere" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(16; R = 1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            @test grid.sin_sg[p, i, j, pos] > 0
        end
    end
end

@testitem "sin_sg: panel symmetry" setup = [FV3Setup] tags = [:fv3] begin
    # On the gnomonic grid, all 6 panels have identical metric properties
    # sin_sg should be the same on every panel at corresponding positions
    grid = CubedSphereGrid(8; R = 1.0)
    for pos in 1:9
        for i in 1:grid.Nc, j in 1:grid.Nc
            vals = [grid.sin_sg[p, i, j, pos] for p in 1:6]
            @test maximum(vals) - minimum(vals) < 1.0e-12
        end
    end
end

@testitem "sin_sg: center matches gnomonic_metric" setup = [FV3Setup] tags = [:fv3] begin
    # Position 5 (cell center) should match the directly computed angle
    grid = CubedSphereGrid(12; R = 6.371e6)
    for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
        J, gxx, gee, gxe = EarthSciDiscretizations.gnomonic_metric(
            grid.ξ_centers[i], grid.η_centers[j], grid.R
        )
        sin_expected = J / sqrt(gxx * gee)
        cos_expected = gxe / sqrt(gxx * gee)
        @test isapprox(grid.sin_sg[p, i, j, 5], sin_expected; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 5], cos_expected; rtol = 1.0e-12)
    end
end

@testitem "sin_sg: shared edge consistency" setup = [FV3Setup] tags = [:fv3] begin
    # The east mid-edge of cell (i,j) is the same point as the west mid-edge of cell (i+1,j)
    # sin_sg[p,i,j,3] (east) should equal sin_sg[p,i+1,j,1] (west) for interior cells
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc
    for p in 1:6, i in 1:(Nc - 1), j in 1:Nc
        @test isapprox(grid.sin_sg[p, i, j, 3], grid.sin_sg[p, i + 1, j, 1]; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 3], grid.cos_sg[p, i + 1, j, 1]; rtol = 1.0e-12)
    end
    # Similarly, north of (i,j) = south of (i,j+1)
    for p in 1:6, i in 1:Nc, j in 1:(Nc - 1)
        @test isapprox(grid.sin_sg[p, i, j, 4], grid.sin_sg[p, i, j + 1, 2]; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 4], grid.cos_sg[p, i, j + 1, 2]; rtol = 1.0e-12)
    end
end

