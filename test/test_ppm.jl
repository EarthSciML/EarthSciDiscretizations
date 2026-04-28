@testsnippet PPMSetup begin
    using Test
    using EarthSciDiscretizations
end

# Loop-form `ppm_reconstruction` was retired in dsc-o05; the 4th-order
# edge interpolation + CW84 limiter stencil is exercised against the
# rule-engine path by `tests/fixtures/ppm_reconstruction/` (Layer-B
# convergence) and per-binding rule conformance suites. The CW84 limiter
# helper itself remains in `src/operators/reconstruction.jl` because
# `flux_1d.jl` and `vertical_remap.jl` (Bucket B, schema-gated) still
# call it; the test below pins its cell-average-preserving contract.

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
