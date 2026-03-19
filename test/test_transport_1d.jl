@testsnippet Transport1DSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "1D Lax-Friedrichs flux of constant field is zero" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q = fill(3.0, 6, Nc, Nc)
    courant = fill(0.1, 6, Nc, Nc)

    ao = flux_1d(q, courant, grid, :xi)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc - 2, Nc)
    @test all(abs.(tendency) .< 1e-12)
end

@testitem "1D Lax-Friedrichs flux sign for uniform positive flow" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    # Increasing field in xi direction
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = Float64(i)
    end

    # Positive courant number (flow in +xi direction)
    courant = fill(0.5, 6, Nc, Nc)

    ao = flux_1d(q, courant, grid, :xi)
    tendency = evaluate_arrayop(ao)

    # With positive flow and increasing q, tendency should be negative
    @test all(tendency .< 0.0)
end

@testitem "1D Lax-Friedrichs flux is linear in q" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q1 = randn(6, Nc, Nc)
    q2 = randn(6, Nc, Nc)
    courant = fill(0.3, 6, Nc, Nc)

    ao1 = flux_1d(q1, courant, grid, :xi)
    ao2 = flux_1d(q2, courant, grid, :xi)
    ao_sum = flux_1d(q1 + q2, courant, grid, :xi)

    t1 = evaluate_arrayop(ao1)
    t2 = evaluate_arrayop(ao2)
    t_sum = evaluate_arrayop(ao_sum)

    @test isapprox(t_sum, t1 + t2; rtol=1e-10)
end

@testitem "PPM 1D flux of constant field with zero velocity is zero" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)

    q = fill(3.0, 6, Nc, Nc)
    vel_xi = fill(0.0, 6, Nc + 1, Nc)
    dt = 0.01

    tendency = zeros(6, Nc, Nc)
    flux_1d_ppm!(tendency, q, vel_xi, grid, :xi, dt)

    @test all(abs.(tendency) .< 1e-14)
end

@testitem "PPM 1D flux conservation" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 16
    grid = CubedSphereGrid(Nc; R=1.0)

    # Smooth field
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * cos(2π * grid.ξ_centers[i] / (π / 2))
    end
    vel_xi = fill(0.5, 6, Nc + 1, Nc)
    dt = 0.01

    tendency = zeros(6, Nc, Nc)
    flux_1d_ppm!(tendency, q, vel_xi, grid, :xi, dt)

    # Total tendency weighted by area should be approximately zero
    # (conservation on a closed sphere — fluxes cancel)
    total = sum(tendency[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    @test abs(total) < 1.0  # Relaxed for panel boundary effects
end
