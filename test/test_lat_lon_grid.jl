@testsnippet LatLonSetup begin
    using Test
    using EarthSciDiscretizations
    const grids = EarthSciDiscretizations.grids
end

# ------------------------------------------------------------------ API --

@testitem "lat_lon: exported from grids namespace" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test isdefined(EarthSciDiscretizations.grids, :lat_lon)
    @test isdefined(EarthSciDiscretizations, :LatLonGrid)
end

@testitem "lat_lon: defaults match GRIDS_API §2.3" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 8, nlat = 4)
    @test g isa LatLonGrid{Float64}
    @test g.variant === :regular
    @test g.R == 6.371e6
    @test g.dtype == "float64"
    @test g.ghosts == 0
    @test g.pole_policy === :none
    @test EarthSciDiscretizations.family(g) == "lat_lon"
    @test nlon_uniform(g) == 8
    @test g.nlat == 4
    @test EarthSciDiscretizations.n_cells(g) == 32
    @test g.lon_start ≈ -pi
end

@testitem "lat_lon: missing nlon / nlat raise ArgumentError" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(; nlat = 4)
    @test_throws ArgumentError grids.lat_lon(; nlon = 4)
    @test_throws ArgumentError grids.lat_lon()
end

@testitem "lat_lon: invalid parameter values rejected" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws DomainError grids.lat_lon(; nlon = 0, nlat = 4)
    @test_throws DomainError grids.lat_lon(; nlon = 4, nlat = 0)
    @test_throws DomainError grids.lat_lon(; nlon = 4, nlat = 4, R = 0)
    @test_throws DomainError grids.lat_lon(; nlon = 4, nlat = 4, R = -1.0)
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 4, ghosts = -1)
end

@testitem "lat_lon: invalid variant / pole_policy rejected" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 4, variant = :bogus)
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 4, pole_policy = :spin)
end

@testitem "lat_lon: non-:none pole policies declared but not implemented" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 4, pole_policy = :average)
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 4, pole_policy = :fold)
end

@testitem "lat_lon: regular variant rejects nlon_per_row" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(; nlon = 4, nlat = 2, nlon_per_row = [4, 4])
end

@testitem "lat_lon: reduced_gaussian rejects nlon" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(;
        variant = :reduced_gaussian, nlon = 4, nlon_per_row = [4, 8]
    )
end

@testitem "lat_lon: reduced_gaussian requires nlon_per_row" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(; variant = :reduced_gaussian, nlat = 4)
end

@testitem "lat_lon: reduced_gaussian nlat / nlon_per_row length mismatch" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws ArgumentError grids.lat_lon(;
        variant = :reduced_gaussian, nlat = 3, nlon_per_row = [4, 8]
    )
end

@testitem "lat_lon: reduced_gaussian rejects zero-width row" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws DomainError grids.lat_lon(;
        variant = :reduced_gaussian, nlon_per_row = [4, 0, 4]
    )
end

@testitem "lat_lon: invalid lat_edges rejected" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    # non-monotonic
    @test_throws DomainError grids.lat_lon(;
        nlon = 2, nlat = 2, lat_edges = [0.0, 0.0, 0.5]
    )
    # outside [-pi/2, pi/2]
    @test_throws DomainError grids.lat_lon(;
        nlon = 2, nlat = 1, lat_edges = [-2.0, 2.0]
    )
end

@testitem "lat_lon: lat_centers outside enclosing edges rejected" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    @test_throws DomainError grids.lat_lon(;
        nlon = 2, nlat = 2,
        lat_edges = [-0.5, 0.0, 0.5],
        lat_centers = [1.0, 0.25],
    )
end

# -------------------------------------------------------------- topology --

@testitem "lat_lon: default lat_edges span pole to pole" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 6)
    @test length(g.lat_edges) == 7
    @test g.lat_edges[1] ≈ -pi / 2
    @test g.lat_edges[end] ≈ pi / 2
end

@testitem "lat_lon: interior neighbours are local" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 8, nlat = 6)
    ns = neighbors(g, 3, 4)
    @test ns[:W] == (3, 3)
    @test ns[:E] == (3, 5)
    @test ns[:S] == (2, 4)
    @test ns[:N] == (4, 4)
end

@testitem "lat_lon: longitude wraps periodically" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 8, nlat = 4)
    @test neighbors(g, 3, 1)[:W] == (3, 8)
    @test neighbors(g, 3, 8)[:E] == (3, 1)
end

@testitem "lat_lon: poles have no N/S neighbour under :none" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 8, nlat = 4)
    ns_s = neighbors(g, 1, 4)
    @test ns_s[:S] === nothing
    @test ns_s[:N] == (2, 4)
    ns_n = neighbors(g, 4, 4)
    @test ns_n[:N] === nothing
    @test ns_n[:S] == (3, 4)
end

@testitem "lat_lon: all neighbours stay in range" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 5, nlat = 4)
    for j in 1:g.nlat
        for i in 1:nlon(g, j)
            for cell in values(neighbors(g, j, i))
                cell === nothing && continue
                jj, ii = cell
                @test 1 ≤ jj ≤ g.nlat
                @test 1 ≤ ii ≤ nlon(g, jj)
            end
        end
    end
end

@testitem "lat_lon: reduced_gaussian neighbour rounding" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    # Row 2 has width 8; row 1 has width 4. Cell (2, 6) center frac =
    # (6 - 0.5) / 8 = 0.6875 → floor(0.6875 * 4) = 2 → 1-based col 3.
    g = grids.lat_lon(; variant = :reduced_gaussian, nlon_per_row = [4, 8, 4])
    @test neighbors(g, 2, 6)[:S] == (1, 3)
    @test neighbors(g, 2, 1)[:S] == (1, 1)
end

@testitem "lat_lon: reduced_gaussian cell counts and row offsets" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; variant = :reduced_gaussian, nlon_per_row = [4, 8, 4])
    @test EarthSciDiscretizations.n_cells(g) == 16
    @test g.nlat == 3
    @test row_offset(g, 1) == 0
    @test row_offset(g, 2) == 4
    @test row_offset(g, 3) == 12
    @test nlon_uniform(g) === nothing
end

@testitem "lat_lon: out-of-range cell index errors" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 4)
    @test_throws ArgumentError neighbors(g, 5, 1)
    @test_throws ArgumentError neighbors(g, 1, 5)
    @test_throws ArgumentError cell_centers(g, 0, 1)
end

# --------------------------------------------------------------- centers --

@testitem "lat_lon: cell_center lat matches row center" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    # nlat=2 → edges at [-pi/2, 0, pi/2]; row 1 center = -pi/4, row 2 = pi/4.
    g = grids.lat_lon(; nlon = 4, nlat = 2)
    _, lat1 = cell_centers(g, 1, 1)
    _, lat2 = cell_centers(g, 2, 1)
    @test lat1 ≈ -pi / 4
    @test lat2 ≈ pi / 4
end

@testitem "lat_lon: bulk cell_centers matches per-cell" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 6, nlat = 4)
    bulk = cell_centers(g)
    for j in (1, 2, 4), i in (1, 3, 6)
        s_lon, s_lat = cell_centers(g, j, i)
        idx = row_offset(g, j) + i
        @test bulk.lon[idx] ≈ s_lon
        @test bulk.lat[idx] ≈ s_lat
    end
end

@testitem "lat_lon: cell_centers lon/lat lie in expected ranges" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 12, nlat = 6)
    bulk = cell_centers(g)
    @test all(bulk.lon .≥ -pi - 1.0e-12)
    @test all(bulk.lon .≤ pi + 1.0e-12)
    @test all(bulk.lat .≥ -pi / 2 - 1.0e-12)
    @test all(bulk.lat .≤ pi / 2 + 1.0e-12)
end

@testitem "lat_lon: lon_edges / lon_centers have correct lengths" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 5, nlat = 3)
    @test length(lon_edges(g, 1)) == 6
    @test length(lon_centers(g, 1)) == 5
end

# ------------------------------------------------------------------ area --

@testitem "lat_lon: total area = 4π R²" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    r = 6.371e6
    g = grids.lat_lon(; nlon = 72, nlat = 36, R = r)
    tot = sum(cell_area(g, j, i) for j in 1:g.nlat for i in 1:nlon(g, j))
    expected = 4.0 * pi * r * r
    @test abs(tot - expected) / expected < 1.0e-12
end

@testitem "lat_lon: total unit-sphere area = 4π" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 10, nlat = 8, R = 1.0)
    tot = sum(cell_area(g, j, i) for j in 1:g.nlat for i in 1:nlon(g, j))
    @test abs(tot - 4.0 * pi) < 1.0e-12
end

@testitem "lat_lon: reduced_gaussian total area = 4π R²" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(;
        variant = :reduced_gaussian,
        nlon_per_row = [4, 8, 12, 12, 8, 4],
        R = 1.0,
    )
    tot = sum(cell_area(g, j, i) for j in 1:g.nlat for i in 1:nlon(g, j))
    @test abs(tot - 4.0 * pi) < 1.0e-12
end

@testitem "lat_lon: area positive everywhere" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 6, nlat = 4)
    for j in 1:g.nlat, i in 1:nlon(g, j)
        @test cell_area(g, j, i) > 0.0
    end
end

@testitem "lat_lon: area via metric_eval matches cell_area" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 6, nlat = 4)
    for j in 1:g.nlat, i in 1:nlon(g, j)
        @test metric_eval(g, :area, j, i) ≈ cell_area(g, j, i)
    end
end

# ---------------------------------------------------------------- metric --

@testitem "lat_lon: metric at equator row ≈ R²" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    r = 3.0
    g = grids.lat_lon(; nlon = 4, nlat = 400, R = r)
    j_eq = 201   # Row nearest the equator (lat_center ≈ 0.00393)
    g_ll = metric_eval(g, :g_lonlon, j_eq, 1)
    g_pp = metric_eval(g, :g_latlat, j_eq, 1)
    @test abs(g_ll - r * r) / (r * r) < 1.0e-4
    @test g_pp ≈ r * r
end

@testitem "lat_lon: g_lonlat = 0 everywhere" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 4)
    for j in 1:g.nlat
        @test metric_eval(g, :g_lonlat, j, 1) == 0.0
        @test metric_eval(g, :ginv_lonlat, j, 1) == 0.0
    end
end

@testitem "lat_lon: inverse metric is inverse" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 6, nlat = 6, R = 1.0)
    for j in (1, 3, 6)
        g_ll = metric_eval(g, :g_lonlon, j, 1)
        g_pp = metric_eval(g, :g_latlat, j, 1)
        i_ll = metric_eval(g, :ginv_lonlon, j, 1)
        i_pp = metric_eval(g, :ginv_latlat, j, 1)
        @test g_ll * i_ll ≈ 1.0
        @test g_pp * i_pp ≈ 1.0
    end
end

@testitem "lat_lon: metric_eval rejects bad index / name" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 4)
    @test_throws ArgumentError metric_eval(g, :J, 5, 1)
    @test_throws ArgumentError metric_eval(g, :J, 1, 5)
    @test_throws ArgumentError metric_eval(g, :not_a_metric, 1, 1)
end

@testitem "lat_lon: Jacobian positive everywhere" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 6, nlat = 6)
    for j in 1:g.nlat, i in 1:nlon(g, j)
        @test metric_eval(g, :J, j, i) > 0.0
    end
end

# ------------------------------------------------------------------ esm --

@testitem "lat_lon: to_esm (regular) is declarative" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 72, nlat = 36, R = 6.371e6, ghosts = 2)
    doc = to_esm(g)
    @test doc["family"] == "lat_lon"
    @test doc["version"] == "1.0.0"
    @test doc["dtype"] == "float64"
    @test doc["topology"] == "rectilinear"
    @test doc["variant"] == "regular"
    @test doc["generator"] == "lat_lon_regular"
    @test doc["params"]["nlon"] == 72
    @test doc["params"]["nlat"] == 36
    @test doc["params"]["R"] == 6.371e6
    @test doc["params"]["ghosts"] == 2
    @test doc["params"]["pole_policy"] == "none"
    @test haskey(doc, "provenance")

    # Mayor 2026-04-20 correction: no inline geometry blobs.
    for forbidden in ("cells", "edges", "vertices", "coordinates", "lon_array")
        @test !haskey(doc, forbidden)
    end
end

@testitem "lat_lon: to_esm (reduced_gaussian) carries schedule" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(;
        variant = :reduced_gaussian, nlon_per_row = [4, 8, 4], R = 1.0
    )
    doc = to_esm(g)
    @test doc["variant"] == "reduced_gaussian"
    @test doc["generator"] == "lat_lon_reduced_gaussian"
    @test doc["params"]["nlat"] == 3
    @test doc["params"]["nlon_per_row"] == [4, 8, 4]
    @test doc["params"]["lat_edges"] isa AbstractVector
end

@testitem "lat_lon: provenance identifies julia binding" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 4)
    prov = to_esm(g)["provenance"]
    @test prov["binding"] == "julia"
    @test prov["source"] == "EarthSciDiscretizations.grids.lat_lon"
    @test prov["generator"] == "lat_lon_regular"
end

@testitem "lat_lon: dtype Float32 propagates to .esm" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    g = grids.lat_lon(; nlon = 4, nlat = 4, dtype = Float32)
    @test g isa LatLonGrid{Float32}
    @test g.dtype == "float32"
    @test to_esm(g)["dtype"] == "float32"
end

# --- Determinism: same inputs → same .esm payload ---

@testitem "lat_lon: determinism — repeated calls produce equal payloads" setup = [LatLonSetup] tags = [:grid, :lat_lon] begin
    opts = (nlon = 8, nlat = 8, R = 6.371e6, ghosts = 1)
    a = to_esm(grids.lat_lon(; opts...))
    b = to_esm(grids.lat_lon(; opts...))
    @test a == b
end
