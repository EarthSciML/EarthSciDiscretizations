@testsnippet PipelineSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "FVCubedSphere construction" setup=[PipelineSetup] tags=[:pipeline] begin
    disc = FVCubedSphere(48)
    @test disc.Nc == 48
    @test disc.Nk == 0
    @test disc.Ng == 3
    @test disc.transport_scheme == :upwind

    disc2 = FVCubedSphere(24; Nk=30, R=1.0, transport=:ppm)
    @test disc2.Nc == 24
    @test disc2.Nk == 30
    @test disc2.R == 1.0
    @test disc2.transport_scheme == :ppm
end

@testitem "Dimension identification" setup=[PipelineSetup] tags=[:pipeline] begin
    @test identify_dimension(:t) == :t
    @test identify_dimension(:lon) == :xi
    @test identify_dimension(:lat) == :eta
    @test identify_dimension(:x) == :xi
    @test identify_dimension(:y) == :eta
    @test identify_dimension(:z) == :vertical
end

@testitem "Initial condition projection" setup=[PipelineSetup] tags=[:pipeline] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R=1.0)

    q = project_initial_condition(grid) do lon, lat
        cos(lat)
    end

    @test size(q) == (6, Nc, Nc)
    @test all(-1.0 .<= q .<= 1.0)
    center = Nc ÷ 2
    @test q[1, center, center] > 0.9
end

@testitem "End-to-end discretize pipeline" setup=[PipelineSetup] tags=[:pipeline] begin
    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using Symbolics
    using DomainSets
    using OrdinaryDiffEqDefault
    using SciMLBase

    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    # Diffusion equation
    eq = [D(u(t, lon, lat)) ~ 0.01 * (Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat))))]
    bcs = [u(0, lon, lat) ~ cos(lat) * cos(lon)]
    domains = [t ∈ Interval(0.0, 0.1),
               lon ∈ Interval(-π, π),
               lat ∈ Interval(-π/2, π/2)]
    @named sys = PDESystem(eq, bcs, domains, [t, lon, lat], [u(t, lon, lat)])

    disc = FVCubedSphere(4; R=1.0)
    prob = SciMLBase.discretize(sys, disc)

    @test prob isa ODEProblem
    @test length(prob.u0) == 6 * 4 * 4  # 6 panels × 4 × 4

    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Initial condition should have max ≈ 1 (cos(0)*cos(0))
    @test maximum(prob.u0) > 0.9
end
