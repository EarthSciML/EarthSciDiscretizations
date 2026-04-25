@testsnippet DuoTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractUnstructuredGrid
    const grids = EarthSciDiscretizations.grids
end

@testitem "duo trait: Tier-C bulk arrays at level 0" setup = [DuoTraitSetup] tags = [:grid, :duo, :trait] begin
    g = grids.duo(; loader = (path = "builtin://icosahedral/0", reader = "auto", check = "strict"), R = 1.0)
    @test g isa AbstractUnstructuredGrid
    @test n_cells(g) == 20    # bare icosahedron
    @test n_dims(g) == 1
    @test axis_names(g) == (:cell,)
    @test length(cell_centers(g, :lon)) == 20
    @test length(cell_centers(g, :lat)) == 20
    @test length(cell_widths(g, :cell)) == 20
    @test cell_volume(g) === g.area
    @test sum(cell_volume(g)) ≈ 4π atol = 1.0e-10
end

@testitem "duo trait: triangular slot neighbours" setup = [DuoTraitSetup] tags = [:grid, :duo, :trait] begin
    g = grids.duo(; loader = (path = "builtin://icosahedral/0",), R = 1.0)
    nb1 = neighbor_indices(g, :cell, 1)
    nb3 = neighbor_indices(g, :cell, 3)
    @test length(nb1) == 20 && length(nb3) == 20
    # Closed icosahedral mesh → all neighbour slots populated, no zeros.
    @test all(nb1 .> 0)
    @test all(nb3 .> 0)
    @test_throws ArgumentError neighbor_indices(g, :cell, 4)
    @test_throws ArgumentError neighbor_indices(g, :cell, 0)
    @test all(.!boundary_mask(g, :cell, :lower))
end

@testitem "duo trait: Tier-U adjacency table" setup = [DuoTraitSetup] tags = [:grid, :duo, :trait] begin
    g = grids.duo(; loader = (path = "builtin://icosahedral/0",), R = 1.0)
    tbl = cell_neighbor_table(g)
    @test size(tbl) == (3, 20)
    val = cell_valence(g)
    @test val == fill(3, 20)
end
