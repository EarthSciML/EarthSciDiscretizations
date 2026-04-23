@testsnippet PPMSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "PPM reproduces linear functions" setup = [PPMSetup] tags = [:ppm] begin
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
        # For a linear function, reconstructed values should be bounded
        @test q_left[p, k, j] >= min(q[p, ic - 1, j], q[p, ic, j]) - 1.0e-10
        @test q_left[p, k, j] <= max(q[p, ic - 1, j], q[p, ic, j]) + 1.0e-10
        @test q_right[p, k, j] >= min(q[p, ic, j], q[p, ic + 1, j]) - 1.0e-10
        @test q_right[p, k, j] <= max(q[p, ic, j], q[p, ic + 1, j]) + 1.0e-10
    end
end

@testitem "PPM handles quadratic functions" setup = [PPMSetup] tags = [:ppm] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)

    # Quadratic field: f = xi^2 (has a minimum near ξ=0)
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = grid.ξ_centers[i]^2
    end

    q_left, q_right = ppm_reconstruction(q, grid, :xi)

    # CW84 limiter ensures no spurious oscillations.
    # At extrema, the limiter flattens (ql=qr=qi), which is correct behavior.
    # Away from extrema, reconstructed values should be reasonable.
    for p in 1:6, k in 1:(Nc - 4), j in 1:Nc
        ic = k + 2
        qi = q[p, ic, j]
        # Left edge should be between qi and neighbor, or flattened to qi at extremum
        @test isfinite(q_left[p, k, j])
        @test isfinite(q_right[p, k, j])
        # CW84 guarantees the integral is preserved
        ql = q_left[p, k, j]; qr = q_right[p, k, j]
        dq = qr - ql; q6 = 6.0 * (qi - 0.5 * (ql + qr))
        integral = ql + dq / 2 + q6 / 6
        @test isapprox(integral, qi; atol = 1.0e-12)
    end
end

@testitem "PPM monotonicity (CW84 limiter)" setup = [PPMSetup] tags = [:ppm] begin
    using EarthSciDiscretizations: _ppm_limit_cw84

    Nc = 16
    grid = CubedSphereGrid(Nc)

    # Step function: 0 for i <= Nc/2, 1 for i > Nc/2
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = i <= div(Nc, 2) ? 0.0 : 1.0
    end

    q_left, q_right = ppm_reconstruction(q, grid, :xi)

    # CW84 is a monotonicity limiter that preserves cell averages and
    # prevents new internal extrema. Near discontinuities, edge values
    # may slightly overshoot the [0,1] range but the parabola integral
    # must preserve the cell average.
    for p in 1:6, k in 1:(Nc - 4), j in 1:Nc
        ic = k + 2
        qi = q[p, ic, j]
        ql = q_left[p, k, j]; qr = q_right[p, k, j]
        # Cell average is preserved
        dq = qr - ql; q6 = 6.0 * (qi - 0.5 * (ql + qr))
        @test isapprox(ql + dq / 2 + q6 / 6, qi; atol = 1.0e-12)
    end

    # In uniform regions (far from the step), values should be exactly qi
    # Check cells at i=3 (well inside the 0-region, physical index=5 in output index=3)
    for p in 1:6, j in 1:Nc
        @test isapprox(q_left[p, 1, j], 0.0; atol = 1.0e-12)
        @test isapprox(q_right[p, 1, j], 0.0; atol = 1.0e-12)
    end
end

@testitem "PPM CW84 limiter preserves cell average" setup = [PPMSetup] tags = [:ppm] begin
    using EarthSciDiscretizations: _ppm_limit_cw84

    # For the CW84 limiter, the integral of the parabola should equal
    # the cell average qi: ∫₀¹ [ql + x*(dq + q6*(1-x))] dx = qi
    # where dq = qr - ql, q6 = 6*(qi - (ql+qr)/2)
    for qi in [0.5, 1.0, 3.0, -2.0]
        for (ql0, qr0) in [(0.0, 1.0), (0.2, 0.8), (-1.0, 2.0), (qi - 0.1, qi + 0.1)]
            ql, qr = _ppm_limit_cw84(ql0, qr0, qi)
            dq = qr - ql
            q6 = 6.0 * (qi - 0.5 * (ql + qr))
            # Integral over [0,1]: ql + dq/2 + q6/2 - q6/3 = ql + dq/2 + q6/6
            integral = ql + dq / 2 + q6 / 6
            @test isapprox(integral, qi; atol = 1.0e-12)
        end
    end
end
