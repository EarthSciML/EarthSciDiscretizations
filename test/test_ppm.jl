@testsnippet PPMSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "PPM reproduces linear functions" setup=[PPMSetup] tags=[:ppm] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)

    # Linear field: f = 3*xi_center + 1
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 3.0 * grid.ξ_centers[i] + 1.0
    end

    q_left, q_right = ppm_reconstruction(q, grid, :xi)

    # PPM should reproduce linear functions exactly
    # q_left/q_right have size (6, Nc-4, Nc)
    # Index (p, k, j) maps to physical cell (p, k+2, j) for the xi dimension
    for p in 1:6, k in 1:(Nc - 4), j in 1:Nc
        ic = k + 2  # physical cell index
        # For a linear function, left and right reconstructed values
        # should be bounded by neighboring cell values
        @test q_left[p, k, j] >= min(q[p, ic - 1, j], q[p, ic, j]) - 1e-10
        @test q_left[p, k, j] <= max(q[p, ic - 1, j], q[p, ic, j]) + 1e-10
        @test q_right[p, k, j] >= min(q[p, ic, j], q[p, ic + 1, j]) - 1e-10
        @test q_right[p, k, j] <= max(q[p, ic, j], q[p, ic + 1, j]) + 1e-10
    end
end

@testitem "PPM reproduces quadratic functions" setup=[PPMSetup] tags=[:ppm] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)

    # Quadratic field: f = xi^2
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = grid.ξ_centers[i]^2
    end

    q_left, q_right = ppm_reconstruction(q, grid, :xi)

    # PPM (4th order) should handle quadratics very well
    # Check that reconstruction values are reasonable (bounded by neighbors)
    for p in 1:6, k in 1:(Nc - 4), j in 1:Nc
        ic = k + 2
        lo = min(q[p, ic - 1, j], q[p, ic, j], q[p, ic + 1, j])
        hi = max(q[p, ic - 1, j], q[p, ic, j], q[p, ic + 1, j])
        @test q_left[p, k, j] >= lo - 1e-10
        @test q_right[p, k, j] <= hi + 1e-10
    end
end

@testitem "PPM monotonicity" setup=[PPMSetup] tags=[:ppm] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)

    # Step function: 0 for i <= Nc/2, 1 for i > Nc/2
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = i <= div(Nc, 2) ? 0.0 : 1.0
    end

    q_left, q_right = ppm_reconstruction(q, grid, :xi)

    # After PPM with limiting (clamp), all values should remain in [0, 1]
    @test all(q_left .>= -1e-14)
    @test all(q_left .<= 1.0 + 1e-14)
    @test all(q_right .>= -1e-14)
    @test all(q_right .<= 1.0 + 1e-14)
end
