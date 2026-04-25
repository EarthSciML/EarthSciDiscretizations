@testsnippet LatLonTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractCurvilinearGrid
    const grids = EarthSciDiscretizations.grids
end

@testitem "latlon trait: regular Tier-C bulk arrays" setup = [LatLonTraitSetup] tags = [:grid, :lat_lon, :trait] begin
    g = grids.lat_lon(; nlon = 8, nlat = 4, R = 1.0)
    @test g isa AbstractCurvilinearGrid
    @test n_cells(g) == 32
    @test n_dims(g) == 2
    @test axis_names(g) == (:lon, :lat)
    cl = cell_centers(g, :lon); cφ = cell_centers(g, :lat)
    @test length(cl) == 32 && length(cφ) == 32
    # Row-major (i fast within row): k=1 is row=1, i=1.
    @test cφ[1] ≈ cφ[8]
    @test cl[1] < cl[8]
    # Cell widths: lon = 2π/n_i, lat = π/4 for nlat=4 with default equal-angle.
    @test cell_widths(g, :lon) ≈ fill(2π / 8, 32)
    @test cell_widths(g, :lat) ≈ fill(π / 4, 32)
    # Volume = R² * dlon * (sin(lat_n) - sin(lat_s)).
    vol = cell_volume(g)
    @test sum(vol) ≈ 4π atol = 1.0e-10
end

@testitem "latlon trait: longitudinal periodic neighbours" setup = [LatLonTraitSetup] tags = [:grid, :lat_lon, :trait] begin
    g = grids.lat_lon(; nlon = 6, nlat = 3, R = 1.0)
    nb_e = neighbor_indices(g, :lon, +1)
    nb_w = neighbor_indices(g, :lon, -1)
    # First row, last cell wraps to first cell of same row.
    @test nb_e[6] == 1
    @test nb_w[1] == 6
    @test nb_e[7] == 8     # row 2, i=1 → i=2
    @test nb_e[12] == 7    # row 2, last cell wraps
end

@testitem "latlon trait: pole sentinels under :none" setup = [LatLonTraitSetup] tags = [:grid, :lat_lon, :trait] begin
    g = grids.lat_lon(; nlon = 4, nlat = 3, R = 1.0)
    # n_cells = 12. Row-1 cells (1..4) have no :S neighbour.
    nb_s = neighbor_indices(g, :lat, -1)
    @test all(nb_s[1:4] .== 0)
    # Row-3 cells (9..12) have no :N neighbour.
    nb_n = neighbor_indices(g, :lat, +1)
    @test all(nb_n[9:12] .== 0)
    # Interior row-2 cells have valid lat neighbours.
    @test all(1 .<= nb_n[5:8] .<= 12)
    @test all(1 .<= nb_s[5:8] .<= 12)
end

@testitem "latlon trait: Tier-M metric tensor" setup = [LatLonTraitSetup] tags = [:grid, :lat_lon, :trait] begin
    g = grids.lat_lon(; nlon = 4, nlat = 2, R = 2.0)
    G = metric_g(g); Ginv = metric_ginv(g)
    @test size(G) == (8, 2, 2)
    # Equator-row latitudes (default equal-angle): lat_centers = (-π/4, π/4)
    cos_lat = cos(π / 4)
    @test G[1, 1, 1] ≈ 4.0 * cos_lat^2 atol = 1.0e-12
    @test G[1, 2, 2] ≈ 4.0
    @test G[1, 1, 2] == 0.0 && G[1, 2, 1] == 0.0
    @test Ginv[1, 1, 1] ≈ 1.0 / G[1, 1, 1] atol = 1.0e-12
    @test Ginv[1, 2, 2] ≈ 1.0 / 4.0
    J = metric_jacobian(g)
    @test J[1] ≈ 4.0 * cos_lat
end

@testitem "latlon trait: reduced-Gaussian ragged layout" setup = [LatLonTraitSetup] tags = [:grid, :lat_lon, :trait] begin
    g = grids.lat_lon(; variant = :reduced_gaussian, nlon_per_row = [4, 8, 4], R = 1.0)
    @test n_cells(g) == 16
    cl = cell_centers(g, :lon)
    cφ = cell_centers(g, :lat)
    @test length(cl) == 16 && length(cφ) == 16
    @test cell_widths(g, :lon) ≈ vcat(fill(2π / 4, 4), fill(2π / 8, 8), fill(2π / 4, 4))
    # neighbor_indices(:lat, +1) at row-1 cell-1 maps to nearest row-2 cell.
    nb_n = neighbor_indices(g, :lat, +1)
    @test 5 <= nb_n[1] <= 12
    # Pole-row cells map to 0 across the boundary.
    @test all(nb_n[13:16] .== 0)
end
