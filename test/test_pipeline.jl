@testsnippet PipelineSetup begin
    using Test
    using EarthSciDiscretizations
    using SymbolicUtils: isarrayop
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

@testitem "Operator registry dispatches correctly" setup=[PipelineSetup] tags=[:pipeline] begin
    Nc = 4
    grid = CubedSphereGrid(Nc; R=1.0)
    disc = FVCubedSphere(Nc; R=1.0)
    registry = build_operator_registry(grid, disc)

    phi = ones(6, Nc, Nc)
    grad = apply_operator(registry, :gradient_xi, phi)
    @test isarrayop(grad)
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
