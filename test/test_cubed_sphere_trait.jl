@testsnippet CubedSphereTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractCurvilinearGrid
end

@testitem "cubed_sphere trait: Tier-C panel layout" setup = [CubedSphereTraitSetup] tags = [:grid, :cubed_sphere, :trait] begin
    Nc = 4
    g = CubedSphereGrid(Nc; R = 1.0)
    @test g isa AbstractCurvilinearGrid
    @test n_cells(g) == 6 * Nc * Nc
    @test n_dims(g) == 2
    @test axis_names(g) == (:xi, :eta)
    @test cell_widths(g, :xi) ≈ fill(g.dξ, n_cells(g))
    @test cell_widths(g, :eta) ≈ fill(g.dη, n_cells(g))
    # cell_centers(:xi) at panel 1, (i=1, j=1) is ξ_centers[1].
    cxi = cell_centers(g, :xi)
    @test cxi[1] ≈ g.ξ_centers[1]
    # n_cells × cell_volume must equal 6 * sum of one panel's areas (i.e. 4πR²).
    @test sum(cell_volume(g)) ≈ 4π atol = 1.0e-6
end

@testitem "cubed_sphere trait: interior neighbour indexing" setup = [CubedSphereTraitSetup] tags = [:grid, :cubed_sphere, :trait] begin
    Nc = 4
    g = CubedSphereGrid(Nc)
    nb_e = neighbor_indices(g, :xi, +1)
    nb_w = neighbor_indices(g, :xi, -1)
    # Panel 1, (i=1, j=1) flat index is 1; +ξ → (i=2, j=1) flat index 2.
    @test nb_e[1] == 2
    @test nb_e[Nc] != 0   # crosses panel edge
    # Panel 1, (i=1, j=1) -ξ crosses to panel 5 (West neighbour).
    @test nb_w[1] != 0
end

@testitem "cubed_sphere trait: cell_centers(:lon/:lat) physical axes" setup = [CubedSphereTraitSetup] tags = [:grid, :cubed_sphere, :trait] begin
    Nc = 4
    g = CubedSphereGrid(Nc; R = 1.0)
    clons = cell_centers(g, :lon)
    clats = cell_centers(g, :lat)
    @test length(clons) == n_cells(g)
    @test length(clats) == n_cells(g)
    # Panel 1, (i=1, j=1) flat index is 1; values must match precomputed lon/lat.
    @test clons[1] ≈ g.lon[1, 1, 1]
    @test clats[1] ≈ g.lat[1, 1, 1]
    # lon in [-π, π], lat in [-π/2, π/2]
    @test all(-π .<= clons .<= π)
    @test all(-π / 2 .<= clats .<= π / 2)
    # Unknown axis still throws
    @test_throws ArgumentError cell_centers(g, :z)
    # Memoization: second call returns identical array object
    @test cell_centers(g, :lon) === cell_centers(g, :lon)
end

@testitem "cubed_sphere trait: Tier-M tensor shapes" setup = [CubedSphereTraitSetup] tags = [:grid, :cubed_sphere, :trait] begin
    Nc = 4
    g = CubedSphereGrid(Nc; R = 1.0)
    G = metric_g(g)
    Ginv = metric_ginv(g)
    @test size(G) == (n_cells(g), 2, 2)
    @test size(Ginv) == (n_cells(g), 2, 2)
    # g_ij and g^{ij} should be approximate inverses (2×2).
    for k in (1, 7, 23)
        det_g = G[k, 1, 1] * G[k, 2, 2] - G[k, 1, 2] * G[k, 2, 1]
        @test det_g > 0
        # g_ij g^{jk} ≈ δ_i^k
        a = G[k, 1, 1] * Ginv[k, 1, 1] + G[k, 1, 2] * Ginv[k, 2, 1]
        b = G[k, 1, 1] * Ginv[k, 1, 2] + G[k, 1, 2] * Ginv[k, 2, 2]
        @test a ≈ 1.0 atol = 1.0e-10
        @test b ≈ 0.0 atol = 1.0e-10
    end
    @test length(metric_jacobian(g)) == n_cells(g)
    @test size(coord_jacobian(g, :lon_lat)) == (n_cells(g), 2, 2)
    @test size(coord_jacobian_second(g, :lon_lat)) == (n_cells(g), 2, 2, 2)
end
