@testsnippet CartesianTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractCurvilinearGrid
    const grids = EarthSciDiscretizations.grids
end

@testitem "cartesian trait: 1D Tier-C bulk shapes" setup = [CartesianTraitSetup] tags = [:grid, :cartesian, :trait] begin
    g = grids.cartesian(; nx = 5, extent = [(0.0, 1.0)])
    @test g isa AbstractCurvilinearGrid
    @test n_cells(g) == 5
    @test n_dims(g) == 1
    @test axis_names(g) == (:x,)
    @test cell_centers(g, :x) ≈ collect(0.1:0.2:0.9)
    @test cell_widths(g, :x) ≈ fill(0.2, 5)
    @test cell_volume(g) ≈ fill(0.2, 5)
    @test neighbor_indices(g, :x, +1) == [2, 3, 4, 5, 0]
    @test neighbor_indices(g, :x, -1) == [0, 1, 2, 3, 4]
    @test boundary_mask(g, :x, :lower) == [true, false, false, false, false]
    @test boundary_mask(g, :x, :upper) == [false, false, false, false, true]
end

@testitem "cartesian trait: 2D Tier-C flat indexing" setup = [CartesianTraitSetup] tags = [:grid, :cartesian, :trait] begin
    g = grids.cartesian(; nx = 3, ny = 2, extent = [(0.0, 3.0), (0.0, 2.0)])
    @test n_cells(g) == 6
    @test n_dims(g) == 2
    @test axis_names(g) == (:x, :y)
    cx = cell_centers(g, :x); cy = cell_centers(g, :y)
    @test length(cx) == 6 && length(cy) == 6
    # Row-major (i fastest): flat index k=1 → (i=1,j=1).
    @test cx[1] ≈ 0.5 && cy[1] ≈ 0.5
    @test cx[3] ≈ 2.5 && cy[3] ≈ 0.5
    @test cx[4] ≈ 0.5 && cy[4] ≈ 1.5
    @test cell_widths(g, :x) ≈ fill(1.0, 6)
    @test cell_widths(g, :y) ≈ fill(1.0, 6)
    @test cell_volume(g) ≈ fill(1.0, 6)
    # Neighbor +x at k=1 (i=1,j=1) → flat index of (i=2,j=1) = 2
    @test neighbor_indices(g, :x, +1)[1] == 2
    @test neighbor_indices(g, :x, +1)[3] == 0      # i=3 → boundary
    @test neighbor_indices(g, :y, +1)[1] == 4      # (i=1,j=2)
    @test neighbor_indices(g, :y, +1)[4] == 0      # j=2 → upper boundary
end

@testitem "cartesian trait: Tier-M identity metric" setup = [CartesianTraitSetup] tags = [:grid, :cartesian, :trait] begin
    g = grids.cartesian(; nx = 4, ny = 3, extent = [(0.0, 4.0), (0.0, 3.0)])
    G = metric_g(g)
    Ginv = metric_ginv(g)
    @test size(G) == (12, 2, 2)
    @test size(Ginv) == (12, 2, 2)
    @test all(G[:, 1, 1] .== 1.0) && all(G[:, 2, 2] .== 1.0)
    @test all(G[:, 1, 2] .== 0.0) && all(G[:, 2, 1] .== 0.0)
    @test Ginv == G
    @test metric_jacobian(g) ≈ fill(1.0, 12)
    @test all(metric_dgij_dxk(g) .== 0.0)
    cj = coord_jacobian(g, :cartesian)
    @test size(cj) == (12, 2, 2)
    @test all(cj[:, 1, 1] .== 1.0) && all(cj[:, 2, 2] .== 1.0)
    @test all(coord_jacobian_second(g, :cartesian) .== 0.0)
end

@testitem "cartesian trait: errors on unknown axis / target" setup = [CartesianTraitSetup] tags = [:grid, :cartesian, :trait] begin
    g = grids.cartesian(; nx = 4, extent = [(0.0, 1.0)])
    @test_throws ArgumentError cell_centers(g, :y)
    @test_throws ArgumentError neighbor_indices(g, :z, +1)
    @test_throws ArgumentError coord_jacobian(g, :lon_lat)
end

@testitem "cartesian trait: bulk-array memo returns the same instance" setup = [CartesianTraitSetup] tags = [:grid, :cartesian, :trait] begin
    g = grids.cartesian(; nx = 6, extent = [(0.0, 1.0)])
    a = cell_centers(g, :x)
    b = cell_centers(g, :x)
    @test a === b   # memoized
    c = cell_volume(g); d = cell_volume(g)
    @test c === d
end
