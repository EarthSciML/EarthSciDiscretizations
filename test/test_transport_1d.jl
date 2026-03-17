@testsnippet Transport1DSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "1D flux of constant field is zero" setup=[Transport1DSetup] tags=[:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q = fill(3.0, 6, Nc, Nc)
    courant = fill(0.1, 6, Nc, Nc)

    ao = flux_1d(q, courant, grid, :xi)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc - 2, Nc)
    @test all(abs.(tendency) .< 1e-12)
end

@testitem "1D flux sign for uniform positive flow" setup=[Transport1DSetup] tags=[:transport] begin
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
    # (advection carries higher values from upstream, reducing local q)
    @test all(tendency .< 0.0)
end

@testitem "1D flux is linear in q" setup=[Transport1DSetup] tags=[:transport] begin
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
