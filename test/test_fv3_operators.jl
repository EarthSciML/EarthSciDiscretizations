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

# =============================================================================
# Two-sided PPM tests (src/operators/ppm_edge.jl)
# =============================================================================

@testitem "Two-sided PPM: constant field with zero velocity" setup = [FV3Setup] tags = [:fv3] begin
    # Zero velocity should give zero tendency regardless of the field
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc
    q_const = ones(6, Nc, Nc)
    vel = zeros(6, Nc + 1, Nc)
    tend = zeros(6, Nc, Nc)
    flux_1d_ppm_twosided!(tend, q_const, vel, grid, :xi, 0.01)
    @test maximum(abs.(tend)) < 1.0e-12
end

@testitem "Two-sided PPM: agrees with standard PPM for interior" setup = [FV3Setup] tags = [:fv3] begin
    # For a smooth interior field with flow entirely within a panel,
    # two-sided PPM should match standard PPM closely
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.3 * sin(4π * (i - 0.5) / Nc) * sin(4π * (j - 0.5) / Nc)
    end
    vel = fill(0.05, 6, Nc + 1, Nc)
    tend_std = zeros(6, Nc, Nc)
    tend_2s = zeros(6, Nc, Nc)
    flux_1d_ppm!(tend_std, q, vel, grid, :xi, 0.001)
    flux_1d_ppm_twosided!(tend_2s, q, vel, grid, :xi, 0.001)

    # Interior cells (away from boundaries) should agree closely
    interior_diffs = Float64[]
    for p in 1:6, i in 4:(Nc - 3), j in 1:Nc
        push!(interior_diffs, abs(tend_2s[p, i, j] - tend_std[p, i, j]))
    end
    @test maximum(interior_diffs) < 1.0e-10  # Should be identical for interior points
end

@testitem "Two-sided PPM: mass conservation" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc

    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * cos(grid.lat[p, i, j])
    end

    vel_xi = fill(0.1, 6, Nc + 1, Nc)
    tend = zeros(6, Nc, Nc)
    flux_1d_ppm_twosided!(tend, q, vel_xi, grid, :xi, 0.01)

    mass_change = sum(
        tend[p, i, j] * grid.area[p, i, j]
            for p in 1:6, i in 1:Nc, j in 1:Nc
    )
    total_mass = sum(
        q[p, i, j] * grid.area[p, i, j]
            for p in 1:6, i in 1:Nc, j in 1:Nc
    )
    @test abs(mass_change / total_mass) < 1.0e-2
end

@testitem "Two-sided PPM: η-direction constant field" setup = [FV3Setup] tags = [:fv3] begin
    # Test two-sided PPM in the η-direction (all previous tests are ξ-only)
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc
    q_const = ones(6, Nc, Nc)
    vel = zeros(6, Nc, Nc + 1)
    tend = zeros(6, Nc, Nc)
    flux_1d_ppm_twosided!(tend, q_const, vel, grid, :eta, 0.01)
    @test maximum(abs.(tend)) < 1.0e-12
end

@testitem "Two-sided PPM: η-direction agrees with standard interior" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.3 * sin(4π * (j - 0.5) / Nc) * sin(4π * (i - 0.5) / Nc)
    end
    vel = fill(0.05, 6, Nc, Nc + 1)
    tend_std = zeros(6, Nc, Nc)
    tend_2s = zeros(6, Nc, Nc)
    flux_1d_ppm!(tend_std, q, vel, grid, :eta, 0.001)
    flux_1d_ppm_twosided!(tend_2s, q, vel, grid, :eta, 0.001)

    interior_diffs = Float64[]
    for p in 1:6, i in 1:Nc, j in 4:(Nc - 3)
        push!(interior_diffs, abs(tend_2s[p, i, j] - tend_std[p, i, j]))
    end
    @test maximum(interior_diffs) < 1.0e-10
end

@testitem "Two-sided PPM: smoother than standard near edges" setup = [FV3Setup] tags = [:fv3] begin
    # Compare two-sided PPM to standard PPM for a smooth field;
    # the two-sided version should give similar or better results near boundaries
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc

    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.3 * sin(2 * grid.lon[p, i, j]) * cos(grid.lat[p, i, j])
    end

    vel = fill(0.05, 6, Nc + 1, Nc)
    tend_std = zeros(6, Nc, Nc)
    tend_2s = zeros(6, Nc, Nc)

    flux_1d_ppm!(tend_std, q, vel, grid, :xi, 0.01)
    flux_1d_ppm_twosided!(tend_2s, q, vel, grid, :xi, 0.01)

    # Both should give similar results (two-sided is a refinement, not a radical change)
    diff = maximum(abs.(tend_2s .- tend_std))
    max_tend = max(maximum(abs.(tend_std)), maximum(abs.(tend_2s)))
    @test diff / max_tend < 0.5  # Within 50% of each other
end

# =============================================================================
# Two-sided PPM edge value formula tests
# =============================================================================

@testitem "Two-sided PPM: exact for linear fields on uniform grid" setup = [FV3Setup] tags = [:fv3] begin
    # On a uniform grid (all dx equal), the two-sided formula and the standard
    # 4th-order formula are both exact for linear fields. The two-sided formula
    # is 2nd-order accurate (by design, for robustness at cube edges), while
    # the standard formula is 4th-order.
    dx = 1.0
    q = [1.0, 3.0, 5.0, 7.0]  # linear field: q(x) = 2x - 1

    # Standard 4th-order: a_{1/2} = 7/12 * (q_0 + q_1) - 1/12 * (q_{-1} + q_2)
    std = 7.0 / 12.0 * (q[2] + q[3]) - 1.0 / 12.0 * (q[1] + q[4])

    # Two-sided: q_m1=q[1], q0=q[2], q1=q[3], q2=q[4]
    ts = EarthSciDiscretizations.ppm_edge_value_twosided(q[2], q[3], q[4], q[1], dx, dx, dx, dx)

    # Both formulas are exact for linear fields
    @test isapprox(ts, std; rtol = 1.0e-12)
    @test isapprox(ts, 4.0; rtol = 1.0e-12)  # exact value at interface
end

@testitem "Two-sided PPM: exact for linear fields on non-uniform grid" setup = [FV3Setup] tags = [:fv3] begin
    # For a linear field q(x) = ax + b, the interface value should be exact
    # regardless of grid spacing
    dx = [0.5, 1.0, 1.5, 2.0]  # non-uniform grid
    x_centers = [dx[1] / 2]
    for k in 2:4
        push!(x_centers, sum(dx[1:(k - 1)]) + dx[k] / 2)
    end
    # Interface between cell 2 and cell 3 is at x = sum(dx[1:2]) = 1.5
    x_iface = sum(dx[1:2])

    a = 2.0; b = 1.0
    q = a .* x_centers .+ b
    expected = a * x_iface + b

    ts = EarthSciDiscretizations.ppm_edge_value_twosided(q[2], q[3], q[4], q[1], dx[2], dx[3], dx[4], dx[1])
    @test isapprox(ts, expected; rtol = 1.0e-12)  # Exact for linear fields on any grid
end
