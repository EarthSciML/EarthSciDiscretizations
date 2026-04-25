@testsnippet ArakawaTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractStaggeredGrid
    const grids = EarthSciDiscretizations.grids
end

@testitem "arakawa trait: Tier-C cartesian-base C-stagger" setup = [ArakawaTraitSetup] tags = [:grid, :arakawa, :trait] begin
    base = EarthSciDiscretizations.CartesianBase(; xlo = 0.0, xhi = 4.0, ylo = 0.0, yhi = 2.0, nx = 4, ny = 2)
    g = grids.arakawa(; base = base, stagger = :C)
    @test g isa AbstractStaggeredGrid
    @test n_cells(g) == 8
    @test n_dims(g) == 2
    @test axis_names(g) == (:x, :y)
    cx = cell_centers(g, :x); cy = cell_centers(g, :y)
    @test length(cx) == 8 && length(cy) == 8
    # Row-major (i fast): k=1 → (i=1,j=1) at center (0.5, 0.5)
    @test cx[1] ≈ 0.5 && cy[1] ≈ 0.5
    @test cx[5] ≈ 0.5 && cy[5] ≈ 1.5    # k=5 → (i=1, j=2)
    @test cell_widths(g, :x) ≈ fill(1.0, 8)
    @test cell_widths(g, :y) ≈ fill(1.0, 8)
    @test cell_volume(g) ≈ fill(1.0, 8)
end

@testitem "arakawa trait: neighbour indices respect interior boundaries" setup = [ArakawaTraitSetup] tags = [:grid, :arakawa, :trait] begin
    base = EarthSciDiscretizations.CartesianBase(; xlo = 0.0, xhi = 3.0, ylo = 0.0, yhi = 2.0, nx = 3, ny = 2)
    g = grids.arakawa(; base = base, stagger = :B)
    nb_x = neighbor_indices(g, :x, +1)
    @test nb_x[1] == 2   # (1,1) → (2,1)
    @test nb_x[3] == 0   # (3,1) → boundary
    nb_y = neighbor_indices(g, :y, +1)
    @test nb_y[1] == 4   # (1,1) → (1,2)
    @test nb_y[4] == 0   # (1,2) → upper-y boundary
    bm = boundary_mask(g, :x, :upper)
    @test bm == [false, false, true, false, false, true]
end
