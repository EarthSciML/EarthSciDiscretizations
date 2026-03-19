@testsnippet Transport2DSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "2D Lax-Friedrichs transport of constant field is zero" setup=[Transport2DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q = fill(7.0, 6, Nc, Nc)
    courant_xi = fill(0.2, 6, Nc, Nc)
    courant_eta = fill(-0.1, 6, Nc, Nc)

    ao = transport_2d(q, courant_xi, courant_eta, grid)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc - 2, Nc - 2)
    @test all(abs.(tendency) .< 1e-12)
end

@testitem "2D Lax-Friedrichs transport is linear in q" setup=[Transport2DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q1 = randn(6, Nc, Nc)
    q2 = randn(6, Nc, Nc)
    courant_xi = fill(0.1, 6, Nc, Nc)
    courant_eta = fill(0.1, 6, Nc, Nc)

    ao1 = transport_2d(q1, courant_xi, courant_eta, grid)
    ao2 = transport_2d(q2, courant_xi, courant_eta, grid)
    ao_sum = transport_2d(q1 + q2, courant_xi, courant_eta, grid)

    t1 = evaluate_arrayop(ao1)
    t2 = evaluate_arrayop(ao2)
    t_sum = evaluate_arrayop(ao_sum)

    @test isapprox(t_sum, t1 + t2; rtol=1e-10)
end

@testitem "Lin-Rood 2D transport of constant field with zero velocity is zero" setup=[Transport2DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)

    q = fill(5.0, 6, Nc, Nc)
    vel_xi = fill(0.0, 6, Nc + 1, Nc)
    vel_eta = fill(0.0, 6, Nc, Nc + 1)
    dt = 0.01

    tendency = zeros(6, Nc, Nc)
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)

    @test all(abs.(tendency) .< 1e-14)
end
