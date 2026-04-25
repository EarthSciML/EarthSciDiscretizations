@testsnippet VerticalTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractVerticalGrid
    const grids = EarthSciDiscretizations.grids
end

@testitem "vertical trait: Tier-C sigma column" setup = [VerticalTraitSetup] tags = [:grid, :vertical, :trait] begin
    g = grids.vertical(; coordinate = :sigma, nz = 4)
    @test g isa AbstractVerticalGrid
    @test n_cells(g) == 4
    @test n_dims(g) == 1
    @test axis_names(g) == (:z,)
    @test length(cell_centers(g, :z)) == 4
    @test cell_widths(g, :z) === g.widths
    @test cell_volume(g) === g.widths
    @test neighbor_indices(g, :z, +1) == [2, 3, 4, 0]
    @test neighbor_indices(g, :z, -1) == [0, 1, 2, 3]
    @test boundary_mask(g, :z, :lower) == [true, false, false, false]
    @test boundary_mask(g, :z, :upper) == [false, false, false, true]
end

@testitem "vertical trait: Tier-V column structure" setup = [VerticalTraitSetup] tags = [:grid, :vertical, :trait] begin
    levels = [0.0, 250.0, 750.0, 2000.0]
    g = grids.vertical(; coordinate = :z, levels = levels)
    @test half_levels(g) == levels
    @test layer_thickness(g) ≈ [250.0, 500.0, 1250.0]
    pc = pressure_coefficients(g)
    @test pc.ak isa Vector{Float64} && isempty(pc.ak)
    @test pc.bk isa Vector{Float64} && isempty(pc.bk)
end

@testitem "vertical trait: hybrid eta pressure_coefficients" setup = [VerticalTraitSetup] tags = [:grid, :vertical, :trait] begin
    nz = 3
    p0 = 1.0e5
    ak = [0.0, 0.0, 0.0, 0.0]
    bk = [1.0, 0.7, 0.3, 0.0]
    g = grids.vertical(; coordinate = :eta, ak = ak, bk = bk, p0 = p0)
    pc = pressure_coefficients(g)
    @test pc.ak == ak
    @test pc.bk == bk
    @test pc.p0 == p0
end

@testitem "vertical trait: rejects unknown axis / side" setup = [VerticalTraitSetup] tags = [:grid, :vertical, :trait] begin
    g = grids.vertical(; coordinate = :sigma, nz = 3)
    @test_throws ArgumentError cell_centers(g, :x)
    @test_throws ArgumentError boundary_mask(g, :z, :middle)
end
