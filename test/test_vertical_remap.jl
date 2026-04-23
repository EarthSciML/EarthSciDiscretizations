@testsnippet VertRemapSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Vertical remap of constant is exact" setup = [VertRemapSetup] tags = [:vertical] begin
    Nk = 10; Nc = 4
    q = fill(3.0, 6, Nc, Nc, Nk)
    dp_old = fill(100.0, 6, Nc, Nc, Nk)
    dp_new = fill(100.0, 6, Nc, Nc, Nk)

    q_new = vertical_remap(q, dp_old, dp_new, Nk)
    @test isapprox(q_new, q; atol = 1.0e-12)
end

@testitem "Vertical remap conservation" setup = [VertRemapSetup] tags = [:vertical] begin
    Nk = 10; Nc = 4
    q = zeros(6, Nc, Nc, Nk)
    dp_old = zeros(6, Nc, Nc, Nk)
    dp_new = zeros(6, Nc, Nc, Nk)

    # Linear profile with varying layer thicknesses
    for p in 1:6, i in 1:Nc, j in 1:Nc, k in 1:Nk
        q[p, i, j, k] = Float64(k) / Nk
        dp_old[p, i, j, k] = 100.0
        dp_new[p, i, j, k] = 100.0 + 10.0 * sin(2π * k / Nk)
    end
    # Normalize dp_new to have same total
    for p in 1:6, i in 1:Nc, j in 1:Nc
        total_old = sum(dp_old[p, i, j, :])
        total_new = sum(dp_new[p, i, j, :])
        dp_new[p, i, j, :] .*= total_old / total_new
    end

    q_new = vertical_remap(q, dp_old, dp_new, Nk)

    # Conservation: sum(q * dp) should be preserved for each column
    for p in 1:6, i in 1:Nc, j in 1:Nc
        mass_old = sum(q[p, i, j, k] * dp_old[p, i, j, k] for k in 1:Nk)
        mass_new = sum(q_new[p, i, j, k] * dp_new[p, i, j, k] for k in 1:Nk)
        @test isapprox(mass_old, mass_new; rtol = 1.0e-10)
    end
end

@testitem "Vertical remap monotonicity" setup = [VertRemapSetup] tags = [:vertical] begin
    Nk = 20; Nc = 2
    q = zeros(6, Nc, Nc, Nk)
    dp_old = fill(100.0, 6, Nc, Nc, Nk)
    dp_new = fill(100.0, 6, Nc, Nc, Nk)

    # Step function
    for p in 1:6, i in 1:Nc, j in 1:Nc, k in 1:Nk
        q[p, i, j, k] = k <= Nk ÷ 2 ? 0.0 : 1.0
    end

    # Slightly different new layers
    for p in 1:6, i in 1:Nc, j in 1:Nc, k in 1:Nk
        dp_new[p, i, j, k] = 100.0 + 5.0 * (-1)^k
    end
    for p in 1:6, i in 1:Nc, j in 1:Nc
        total = sum(dp_new[p, i, j, :])
        dp_new[p, i, j, :] .*= sum(dp_old[p, i, j, :]) / total
    end

    q_new = vertical_remap(q, dp_old, dp_new, Nk)

    # CW84 limiter allows small overshoots near discontinuities (~0.5%)
    # while preserving conservation. Values should be close to [0, 1].
    @test all(q_new .>= -0.01)
    @test all(q_new .<= 1.01)
end
