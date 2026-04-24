@testsnippet MpasSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: grids

    # Synthetic tetrahedral Voronoi mesh: 4 triangular cells at the four
    # tetrahedral directions. Each cell has 3 neighbors (the other three);
    # per-cell area = π R² so total area = 4π R². 6 edges; each edge is a
    # great-circle arc between two cell centers. Mirrors the tetra_mesh in
    # python/tests/test_mpas.py but with Julia 1-based indices and `0`
    # sentinel.
    function _tetra_mesh(; R::Float64 = 1.0)
        verts = [
            (1.0, 1.0, 1.0),
            (1.0, -1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ]
        n_cells = 4
        max_edges = 3
        lon_cell = Vector{Float64}(undef, n_cells)
        lat_cell = Vector{Float64}(undef, n_cells)
        x_cell = Vector{Float64}(undef, n_cells)
        y_cell = Vector{Float64}(undef, n_cells)
        z_cell = Vector{Float64}(undef, n_cells)
        for (c, (x, y, z)) in enumerate(verts)
            n = sqrt(x^2 + y^2 + z^2)
            xh, yh, zh = x / n, y / n, z / n
            lon_cell[c] = atan(yh, xh)
            lat_cell[c] = asin(zh)
            x_cell[c] = R * xh
            y_cell[c] = R * yh
            z_cell[c] = R * zh
        end

        area_cell = fill(π * R^2, n_cells)
        n_edges_on_cell = fill(3, n_cells)

        # 1-based adjacency, 0 sentinel.
        cells_on_cell = zeros(Int, max_edges, n_cells)
        for c in 1:n_cells
            j = 0
            for k in 1:n_cells
                k == c && continue
                j += 1
                cells_on_cell[j, c] = k
            end
        end

        n_edges = 6
        cells_on_edge = zeros(Int, 2, n_edges)
        edge_lookup = Dict{Tuple{Int, Int}, Int}()
        e = 0
        for i in 1:n_cells, j in (i + 1):n_cells
            e += 1
            cells_on_edge[1, e] = i
            cells_on_edge[2, e] = j
            edge_lookup[(i, j)] = e
        end

        edges_on_cell = zeros(Int, max_edges, n_cells)
        for c in 1:n_cells, j in 1:max_edges
            nb = cells_on_cell[j, c]
            key = c < nb ? (c, nb) : (nb, c)
            edges_on_cell[j, c] = edge_lookup[key]
        end

        lon_edge = Vector{Float64}(undef, n_edges)
        lat_edge = Vector{Float64}(undef, n_edges)
        dc_edge = Vector{Float64}(undef, n_edges)
        for ei in 1:n_edges
            c1 = cells_on_edge[1, ei]
            c2 = cells_on_edge[2, ei]
            v1 = (x_cell[c1] / R, y_cell[c1] / R, z_cell[c1] / R)
            v2 = (x_cell[c2] / R, y_cell[c2] / R, z_cell[c2] / R)
            d = clamp(v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3], -1.0, 1.0)
            dc_edge[ei] = R * acos(d)
            mx = 0.5 * (v1[1] + v2[1])
            my = 0.5 * (v1[2] + v2[2])
            mz = 0.5 * (v1[3] + v2[3])
            nn = sqrt(mx^2 + my^2 + mz^2)
            lon_edge[ei] = atan(my / nn, mx / nn)
            lat_edge[ei] = asin(mz / nn)
        end
        dv_edge = copy(dc_edge)

        return mpas_mesh_data(;
            lon_cell = lon_cell,
            lat_cell = lat_cell,
            area_cell = area_cell,
            n_edges_on_cell = n_edges_on_cell,
            cells_on_cell = cells_on_cell,
            edges_on_cell = edges_on_cell,
            lon_edge = lon_edge,
            lat_edge = lat_edge,
            cells_on_edge = cells_on_edge,
            dc_edge = dc_edge,
            dv_edge = dv_edge,
            max_edges = max_edges,
            x_cell = x_cell,
            y_cell = y_cell,
            z_cell = z_cell,
            n_vertices = n_cells,
            R = R,
        )
    end
end

@testitem "MPAS: grids.mpas namespace binding" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh(), R = 1.0)
    @test g isa MpasGrid
    @test family(g) == "mpas"
    @test n_cells(g) == 4
    @test n_edges(g) == 6
    @test max_edges(g) == 3
    @test g.ghosts == 0
    @test g.loader === nothing
end

@testitem "MPAS: missing mesh and loader raises" setup = [MpasSetup] tags = [:grid, :mpas] begin
    @test_throws ArgumentError grids.mpas()
    @test_throws ArgumentError grids.mpas(loader = (path = "x",), reader_fn = nothing)
end

@testitem "MPAS: rejects invalid dtype / ghosts / R" setup = [MpasSetup] tags = [:grid, :mpas] begin
    @test_throws ArgumentError grids.mpas(mesh = _tetra_mesh(), dtype = Float16)
    @test_throws DomainError grids.mpas(mesh = _tetra_mesh(), ghosts = 1)
    @test_throws DomainError grids.mpas(mesh = _tetra_mesh(), R = 0.0)
    @test_throws DomainError grids.mpas(mesh = _tetra_mesh(), R = -1.0)
    @test_throws DomainError grids.mpas(mesh = _tetra_mesh(), R = NaN)
end

@testitem "MPAS: topology — neighbor symmetry" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh())
    for c in 1:n_cells(g)
        for nb in neighbors(g, c)
            @test c in neighbors(g, nb)
        end
    end
end

@testitem "MPAS: total area ≈ 4πR²" setup = [MpasSetup] tags = [:grid, :mpas] begin
    for R in (1.0, 6.371e6)
        g = grids.mpas(mesh = _tetra_mesh(R = R), R = R)
        @test isapprox(total_area(g), 4π * R^2; rtol = 1.0e-12)
    end
end

@testitem "MPAS: edge_length great-circle on unit sphere" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh(R = 1.0), R = 1.0)
    # Tetrahedral arc = acos(-1/3)
    expected = acos(-1 / 3)
    for e in 1:n_edges(g)
        @test isapprox(edge_length(g, e), expected; rtol = 1.0e-12)
        @test isapprox(cell_distance(g, e), expected; rtol = 1.0e-12)
    end
end

@testitem "MPAS: cell_centers accessor" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh())
    lon_lat = cell_centers(g)
    @test length(lon_lat.lon) == n_cells(g)
    @test length(lon_lat.lat) == n_cells(g)
    @test all(-π .≤ lon_lat.lon .≤ π)
    @test all(-π / 2 .≤ lon_lat.lat .≤ π / 2)
    lon1, lat1 = cell_centers(g, 1)
    @test isfinite(lon1) && isfinite(lat1)
end

@testitem "MPAS: metric_eval names" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh(R = 1.0), R = 1.0)
    @test isapprox(metric_eval(g, :area, 1), π; rtol = 1.0e-12)
    @test metric_eval(g, "area", 1) == metric_eval(g, :area, 1)
    @test metric_eval(g, :n_edges_on_cell, 1) == 3.0
    @test isapprox(metric_eval(g, :dv_edge, 1), acos(-1 / 3); rtol = 1.0e-12)
    @test_throws ArgumentError metric_eval(g, :bogus, 1)
end

@testitem "MPAS: bounds checks on accessors" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh())
    @test_throws BoundsError cell_centers(g, 0)
    @test_throws BoundsError cell_centers(g, n_cells(g) + 1)
    @test_throws BoundsError cell_area(g, 0)
    @test_throws BoundsError edge_length(g, 0)
    @test_throws BoundsError edge_length(g, n_edges(g) + 1)
end

@testitem "MPAS: loader coercion + strict validation" setup = [MpasSetup] tags = [:grid, :mpas] begin
    # Missing loader.path
    @test_throws ArgumentError EarthSciDiscretizations._coerce_mpas_loader(
        (reader = "auto",)
    )
    # Invalid reader/check
    @test_throws DomainError EarthSciDiscretizations._coerce_mpas_loader(
        (path = "x.nc", reader = "bogus",)
    )
    @test_throws DomainError EarthSciDiscretizations._coerce_mpas_loader(
        (path = "x.nc", check = "sloppy",)
    )
    # NamedTuple and Dict both accepted
    ldr1 = EarthSciDiscretizations._coerce_mpas_loader(
        (path = "mesh.nc", reader = "mpas_mesh", check = "strict")
    )
    @test ldr1 isa MpasLoader
    ldr2 = EarthSciDiscretizations._coerce_mpas_loader(
        Dict("path" => "mesh.nc", "reader" => "auto", "check" => "lenient")
    )
    @test ldr2.reader == "auto"
    @test ldr2.check == "lenient"
end

@testitem "MPAS: loader + reader_fn path" setup = [MpasSetup] tags = [:grid, :mpas] begin
    mesh = _tetra_mesh()
    reader_fn = path -> begin
        @test path == "builtin://tetra"
        return mesh
    end
    g = grids.mpas(
        loader = (path = "builtin://tetra", reader = "mpas_mesh", check = "strict"),
        reader_fn = reader_fn,
    )
    @test g.loader isa MpasLoader
    @test g.loader.path == "builtin://tetra"
    @test n_cells(g) == 4
end

@testitem "MPAS: to_esm declarative shape (no inline geometry)" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh(), R = 6.371e6)
    esm = to_esm(g)
    @test esm["family"] == "mpas"
    @test esm["topology"] == "unstructured"
    @test esm["version"] == "1.0.0"
    @test esm["schema_version"] == "1.0.0"
    @test esm["n_cells"] == 4
    @test esm["n_edges"] == 6
    @test esm["max_edges"] == 3
    @test esm["options"]["R"] == 6.371e6
    @test esm["options"]["loader"] === nothing
    @test haskey(esm["provenance"], "binding")
    @test esm["provenance"]["family"] == "mpas"
    # No inline geometry arrays — declarative config only per dsc-3nw.
    for k in ("lon_cell", "lat_cell", "area_cell", "cells_on_cell")
        @test !haskey(esm, k)
    end
end

@testitem "MPAS: check_mesh strict enforces reciprocity" setup = [MpasSetup] tags = [:grid, :mpas] begin
    mesh = _tetra_mesh()
    # Break reciprocity: cell 1 → 2, but rewrite cell 2's neighbor list to
    # omit 1.
    coc = copy(mesh.cells_on_cell)
    # Replace any entry in cell 2's neighbors that equals 1 with a different cell.
    for j in 1:size(coc, 1)
        if coc[j, 2] == 1
            coc[j, 2] = 3  # duplicate — still in-bounds but breaks symmetry
        end
    end
    broken = MpasMeshData(
        mesh.n_cells, mesh.n_edges, mesh.n_vertices, mesh.max_edges,
        mesh.lon_cell, mesh.lat_cell, mesh.x_cell, mesh.y_cell, mesh.z_cell,
        mesh.area_cell, mesh.n_edges_on_cell, coc, mesh.edges_on_cell,
        mesh.lon_edge, mesh.lat_edge, mesh.cells_on_edge, mesh.dc_edge,
        mesh.dv_edge,
    )
    @test_throws AssertionError check_mesh(broken, true)
    # Lenient skips reciprocity.
    @test check_mesh(broken, false) === nothing
end

@testitem "MPAS: check_mesh bounds checks catch out-of-range indices" setup = [MpasSetup] tags = [:grid, :mpas] begin
    mesh = _tetra_mesh()
    coc = copy(mesh.cells_on_cell)
    coc[1, 1] = 999
    broken = MpasMeshData(
        mesh.n_cells, mesh.n_edges, mesh.n_vertices, mesh.max_edges,
        mesh.lon_cell, mesh.lat_cell, mesh.x_cell, mesh.y_cell, mesh.z_cell,
        mesh.area_cell, mesh.n_edges_on_cell, coc, mesh.edges_on_cell,
        mesh.lon_edge, mesh.lat_edge, mesh.cells_on_edge, mesh.dc_edge,
        mesh.dv_edge,
    )
    @test_throws AssertionError check_mesh(broken, false)
end

@testitem "MPAS: mesh_data shape errors" setup = [MpasSetup] tags = [:grid, :mpas] begin
    @test_throws DomainError mpas_mesh_data(;
        lon_cell = [0.0], lat_cell = [0.0], area_cell = [1.0],
        n_edges_on_cell = [0],
        cells_on_cell = zeros(Int, 3, 1), edges_on_cell = zeros(Int, 3, 1),
        lon_edge = Float64[], lat_edge = Float64[],
        cells_on_edge = zeros(Int, 2, 0),
        dc_edge = Float64[], dv_edge = Float64[],
        max_edges = 0,
    )
    @test_throws ArgumentError mpas_mesh_data(;
        lon_cell = [0.0, 1.0], lat_cell = [0.0],  # mismatched length
        area_cell = [1.0, 1.0], n_edges_on_cell = [0, 0],
        cells_on_cell = zeros(Int, 3, 2), edges_on_cell = zeros(Int, 3, 2),
        lon_edge = Float64[], lat_edge = Float64[],
        cells_on_edge = zeros(Int, 2, 0),
        dc_edge = Float64[], dv_edge = Float64[],
        max_edges = 3,
    )
end

@testitem "MPAS: dtype Float32" setup = [MpasSetup] tags = [:grid, :mpas] begin
    g = grids.mpas(mesh = _tetra_mesh(), dtype = Float32)
    @test g.dtype == "float32"
    @test g.R isa Float32
    @test to_esm(g)["dtype"] == "float32"
end
