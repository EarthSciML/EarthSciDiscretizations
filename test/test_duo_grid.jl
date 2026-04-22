@testsnippet DuoSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: grids
    _builtin(level) = (path = "builtin://icosahedral/$level",
                       reader = "builtin_icosahedral",
                       check = "strict")
end

@testitem "DUO level 0 topology (bare icosahedron)" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(0), R = 1.0)
    @test n_cells(g) == 20       # icosahedron has 20 faces
    @test n_vertices(g) == 12    # 12 vertices
    @test n_edges(g) == 30       # 30 edges
    # Euler characteristic of a sphere: V - E + F = 2
    @test n_vertices(g) - n_edges(g) + n_cells(g) == 2
    @test family(g) == "duo"
end

@testitem "DUO subdivision counts match 20·4^r / 10·4^r+2" setup=[DuoSetup] tags=[:grid, :duo] begin
    for r in 0:3
        g = grids.duo(loader = _builtin(r), R = 1.0)
        @test n_cells(g)    == 20 * 4^r
        @test n_vertices(g) == 10 * 4^r + 2
        @test n_edges(g)    == 30 * 4^r
        @test n_vertices(g) - n_edges(g) + n_cells(g) == 2
    end
end

@testitem "DUO total area ≈ 4πR²" setup=[DuoSetup] tags=[:grid, :duo] begin
    for r in 0:3
        g = grids.duo(loader = _builtin(r), R = 1.0)
        @test isapprox(total_area(g), 4π; rtol = 1e-10)
    end
    R = 6.371e6
    g = grids.duo(loader = _builtin(2), R = R)
    @test isapprox(total_area(g), 4π * R^2; rtol = 1e-10)
end

@testitem "DUO all cells have positive area" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(2), R = 1.0)
    @test all(g.area .> 0)
end

@testitem "DUO cell_neighbors are symmetric and reference valid cells" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(2), R = 1.0)
    Nc = n_cells(g)
    # Closed mesh: every cell has exactly 3 neighbors (no zeros).
    @test all(g.cell_neighbors .> 0)
    @test all(g.cell_neighbors .<= Nc)
    # Reciprocity: if c → n via some edge, then n → c via some edge.
    for c in 1:Nc
        for k in 1:3
            nb = g.cell_neighbors[k, c]
            @test c in (g.cell_neighbors[1, nb], g.cell_neighbors[2, nb], g.cell_neighbors[3, nb])
        end
    end
end

@testitem "DUO neighbors(c) accessor matches cell_neighbors" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(1), R = 1.0)
    for c in 1:n_cells(g)
        nbrs = neighbors(g, c)
        @test nbrs == (g.cell_neighbors[1, c], g.cell_neighbors[2, c], g.cell_neighbors[3, c])
    end
end

@testitem "DUO cell_centers accessor returns lon,lat on sphere" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(1), R = 1.0)
    lon_lat = cell_centers(g)
    @test length(lon_lat.lon) == n_cells(g)
    @test length(lon_lat.lat) == n_cells(g)
    @test all(-π .≤ lon_lat.lon .≤ π)
    @test all(-π/2 .≤ lon_lat.lat .≤ π/2)
    lon1, lat1 = cell_centers(g, 1)
    @test isfinite(lon1) && isfinite(lat1)
end

@testitem "DUO metric_eval delivers area/lon/lat/x/y/z" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(1), R = 6.371e6)
    for c in (1, 5, n_cells(g))
        @test metric_eval(g, :area, c) == g.area[c]
        @test metric_eval(g, :lon, c) == g.lon[c]
        @test metric_eval(g, :lat, c) == g.lat[c]
        x = metric_eval(g, :x, c); y = metric_eval(g, :y, c); z = metric_eval(g, :z, c)
        @test isapprox(sqrt(x^2 + y^2 + z^2), g.R; rtol = 1e-12)
    end
    @test_throws ArgumentError metric_eval(g, :nonsense, 1)
end

@testitem "DUO to_esm returns §6-style declarative config (no inline geometry)" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(1), R = 6.371e6)
    esm = to_esm(g)
    @test esm["family"] == "duo"
    @test esm["topology"] == "unstructured"
    @test esm["dtype"] == "float64"
    @test esm["n_cells"] == n_cells(g)
    @test esm["n_vertices"] == n_vertices(g)
    @test esm["n_edges"] == n_edges(g)
    @test esm["options"]["level"] == 1
    @test esm["options"]["loader"]["path"] == "builtin://icosahedral/1"
    @test esm["options"]["loader"]["reader"] == "builtin_icosahedral"
    @test haskey(esm, "provenance")
    # Mayor's correction: no inline cells/edges/vertices arrays.
    for forbidden in ("cells", "vertices_xyz", "edges_list", "faces")
        @test !haskey(esm, forbidden)
    end
end

@testitem "DUO Float32 dtype" setup=[DuoSetup] tags=[:grid, :duo] begin
    g = grids.duo(loader = _builtin(1), R = 1.0f0, dtype = Float32)
    @test g.dtype == "float32"
    @test eltype(g.area) === Float32
    @test eltype(g.lon) === Float32
    @test isapprox(total_area(g), 4π; rtol = 1e-4)
end

@testitem "DUO invalid loader rejected" setup=[DuoSetup] tags=[:grid, :duo] begin
    # duo_mesh reader without file spec must raise.
    @test_throws ArgumentError grids.duo(loader = (path = "/tmp/no.duo", reader = "duo_mesh", check = "strict"))
    # Negative level rejected.
    @test_throws DomainError grids.duo(loader = _builtin(-1), R = 1.0)
    # Unrecognized path scheme.
    @test_throws ArgumentError grids.duo(loader = (path = "bogus://x", reader = "unknown_reader", check = "strict"))
    # Invalid dtype.
    @test_throws ArgumentError grids.duo(loader = _builtin(0), dtype = Int)
end

@testitem "DUO DuoLoader accepts NamedTuple/Dict/struct" setup=[DuoSetup] tags=[:grid, :duo] begin
    # NamedTuple
    g1 = grids.duo(loader = _builtin(1), R = 1.0)
    # Dict (string keys)
    g2 = grids.duo(loader = Dict("path" => "builtin://icosahedral/1",
                                 "reader" => "builtin_icosahedral",
                                 "check" => "strict"), R = 1.0)
    # Struct
    ldr = DuoLoader(path = "builtin://icosahedral/1", reader = "builtin_icosahedral")
    g3 = grids.duo(loader = ldr, R = 1.0)
    @test n_cells(g1) == n_cells(g2) == n_cells(g3)
    @test g1.lon ≈ g2.lon
    @test g1.lon ≈ g3.lon
end

@testitem "DUO determinism: two calls produce identical grids" setup=[DuoSetup] tags=[:grid, :duo] begin
    g1 = grids.duo(loader = _builtin(2), R = 1.0)
    g2 = grids.duo(loader = _builtin(2), R = 1.0)
    @test g1.lon == g2.lon
    @test g1.lat == g2.lat
    @test g1.area == g2.area
    @test g1.faces == g2.faces
    @test g1.cell_neighbors == g2.cell_neighbors
end

@testitem "DUO cell centers lie on the sphere of radius R" setup=[DuoSetup] tags=[:grid, :duo] begin
    R = 6.371e6
    g = grids.duo(loader = _builtin(2), R = R)
    for c in 1:n_cells(g)
        x = g.cell_cart[1, c]; y = g.cell_cart[2, c]; z = g.cell_cart[3, c]
        @test isapprox(sqrt(x^2 + y^2 + z^2), R; rtol = 1e-12)
    end
end
