@testsnippet MpasConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: grids
    using JSON

    # This conformance harness lives at tests/conformance/grids/mpas/ and
    # is binding-neutral: fixtures.json carries the full MPAS mesh arrays
    # inline (since NetCDF mesh files are not bundled in this repo). Each
    # binding constructs an MpasMeshData from those arrays, builds a
    # grid, queries accessors at pinned points, and compares to the
    # committed golden.

    const MPAS_HARNESS_DIR =
        joinpath(@__DIR__, "..", "tests", "conformance", "grids", "mpas")

    # Python/Rust/TS bindings use 0-based cell/edge indices with -1 as the
    # "no neighbor" / "no edge" sentinel. Julia uses 1-based indices with
    # 0 as the sentinel. The mapping `julia_idx = python_idx + 1` covers
    # both: -1 → 0 (sentinel), 0..N-1 → 1..N (real index).
    _to_julia_idx(i) = Int(i) + 1

    function _build_mpas_grid_from_fixture(fixture::AbstractDict)
        m = fixture["mesh"]
        opts = fixture["opts"]
        max_edges = Int(m["max_edges"])

        lon_cell = Float64.(m["lon_cell"])
        n_cells_local = length(lon_cell)
        lat_cell = Float64.(m["lat_cell"])
        area_cell = Float64.(m["area_cell"])
        n_edges_on_cell = Int.(m["n_edges_on_cell"])

        # Flat row-major arrays in fixtures.json: outer axis = cell, inner
        # axis = slot. Reshape into Julia's column-major (slot, cell)
        # layout, then shift indices.
        flat_coc = Int.(m["cells_on_cell"])
        cells_on_cell = Matrix{Int}(undef, max_edges, n_cells_local)
        @inbounds for c in 1:n_cells_local, j in 1:max_edges
            cells_on_cell[j, c] =
                _to_julia_idx(flat_coc[(c - 1) * max_edges + j])
        end

        flat_eoc = Int.(m["edges_on_cell"])
        edges_on_cell = Matrix{Int}(undef, max_edges, n_cells_local)
        @inbounds for c in 1:n_cells_local, j in 1:max_edges
            edges_on_cell[j, c] =
                _to_julia_idx(flat_eoc[(c - 1) * max_edges + j])
        end

        lon_edge = Float64.(m["lon_edge"])
        n_edges_local = length(lon_edge)
        lat_edge = Float64.(m["lat_edge"])

        flat_coe = Int.(m["cells_on_edge"])
        cells_on_edge = Matrix{Int}(undef, 2, n_edges_local)
        @inbounds for e in 1:n_edges_local, s in 1:2
            cells_on_edge[s, e] =
                _to_julia_idx(flat_coe[(e - 1) * 2 + s])
        end

        dc_edge = Float64.(m["dc_edge"])
        dv_edge = Float64.(m["dv_edge"])

        x_cell = Float64.(m["x_cell"])
        y_cell = Float64.(m["y_cell"])
        z_cell = Float64.(m["z_cell"])
        n_vertices = Int(m["n_vertices"])

        mesh = mpas_mesh_data(;
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
            n_vertices = n_vertices,
            R = Float64(opts["R"]),
        )

        dtype = String(opts["dtype"]) == "float32" ? Float32 : Float64
        return grids.mpas(;
            mesh = mesh,
            R = Float64(opts["R"]),
            dtype = dtype,
            ghosts = Int(opts["ghosts"]),
        )
    end

    # Schema's relative tolerance applied with the same scaling as the
    # python harness: `|a - b| ≤ tol * max(1, |a|, |b|)`.
    function _close_rel(a::Real, b::Real, tol::Real)
        scale = max(1.0, abs(Float64(a)), abs(Float64(b)))
        return abs(Float64(a) - Float64(b)) <= tol * scale
    end

    const _CELL_METRICS =
        ("lon", "lat", "area", "x", "y", "z", "n_edges_on_cell")
    const _EDGE_METRICS = ("lon_edge", "lat_edge", "dc_edge", "dv_edge")
end

@testitem "mpas cross-language conformance" setup = [MpasConformanceSetup] tags =
    [:conformance, :grid, :mpas] begin

    # Always assert the harness files are present so structural drift is
    # caught even before any per-fixture math runs.
    fixtures_path = joinpath(MPAS_HARNESS_DIR, "fixtures.json")
    @test isfile(fixtures_path)
    spec = JSON.parsefile(fixtures_path)
    @test spec["family"] == "mpas"
    @test haskey(spec, "tolerance")
    @test haskey(spec, "fixtures")

    rel_tol = Float64(spec["tolerance"]["relative"])

    for fixture in spec["fixtures"]
        name = fixture["name"]
        golden_path = joinpath(MPAS_HARNESS_DIR, "golden", "$name.json")
        @test isfile(golden_path)

        grid = _build_mpas_grid_from_fixture(fixture)
        golden = JSON.parsefile(golden_path)

        @test n_cells(grid) == golden["n_cells"]
        @test n_edges(grid) == golden["n_edges"]
        @test max_edges(grid) == golden["max_edges"]

        for (k, c0) in enumerate(fixture["query_cells"])
            c = _to_julia_idx(c0)

            lon, lat = cell_centers(grid, c)
            gc = golden["cell_centers"][k]
            @test _close_rel(lon, gc["lon"], rel_tol)
            @test _close_rel(lat, gc["lat"], rel_tol)

            x, y, z = cell_center_cart(grid, c)
            gx = golden["cell_centers_cart"][k]
            @test _close_rel(x, gx["x"], rel_tol)
            @test _close_rel(y, gx["y"], rel_tol)
            @test _close_rel(z, gx["z"], rel_tol)

            got_nb = neighbors(grid, c)
            exp_nb = [_to_julia_idx(n) for n in golden["neighbors"][k]]
            @test got_nb == exp_nb

            @test _close_rel(
                cell_area(grid, c), Float64(golden["cell_area"][k]), rel_tol
            )

            for metric_name in _CELL_METRICS
                got = metric_eval(grid, metric_name, c)
                expected = Float64(golden["cell_metrics"][metric_name][k])
                @test _close_rel(got, expected, rel_tol)
            end
        end

        for (k, e0) in enumerate(fixture["query_edges"])
            e = _to_julia_idx(e0)

            @test _close_rel(
                edge_length(grid, e),
                Float64(golden["edge_length"][k]),
                rel_tol,
            )
            @test _close_rel(
                cell_distance(grid, e),
                Float64(golden["cell_distance"][k]),
                rel_tol,
            )
            for metric_name in _EDGE_METRICS
                got = metric_eval(grid, metric_name, e)
                expected = Float64(golden["edge_metrics"][metric_name][k])
                @test _close_rel(got, expected, rel_tol)
            end
        end
    end
end
