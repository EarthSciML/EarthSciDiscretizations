@testsnippet DuoConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON

    # Julia's DuoGrid is 1-indexed and uses `0` for "no neighbor"; the
    # binding-neutral golden stores 0-based cell indices and `-1` for
    # boundary. Convert at the harness boundary so all four bindings
    # compare byte-for-byte on the same numbers.
    _neighbors_0idx(g, c1::Int) = [
        n == 0 ? -1 : n - 1 for n in EarthSciDiscretizations.neighbors(g, c1)
    ]
end

@testitem "DUO cross-language conformance" setup = [DuoConformanceSetup] tags = [:conformance, :grid, :duo] begin
    HARNESS_DIR = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "duo")
    FIXTURES = JSON.parsefile(joinpath(HARNESS_DIR, "fixtures.json"))
    REL_TOL = Float64(FIXTURES["tolerance"]["relative"])
    METRIC_NAMES = (:area, :lon, :lat, :x, :y, :z)

    function close_rel(a::Real, b::Real, tol::Real)
        scale = max(1.0, abs(a), abs(b))
        return abs(a - b) <= tol * scale
    end

    for fixture in FIXTURES["fixtures"]
        name = fixture["name"]
        opts = fixture["opts"]
        level = Int(opts["level"])
        R = float(opts["R"])
        ghosts = Int(opts["ghosts"])
        golden = JSON.parsefile(joinpath(HARNESS_DIR, "golden", "$name.json"))

        @testset "fixture=$name" begin
            loader = DuoLoader(path = "builtin://icosahedral/$level")
            grid = build_duo_grid(; loader = loader, R = R, ghosts = ghosts)

            @test n_cells(grid) == golden["n_cells"]
            @test n_vertices(grid) == golden["n_vertices"]
            @test n_edges(grid) == golden["n_edges"]

            qps = Int.(fixture["query_points"])
            @test length(qps) == length(golden["cell_centers"])

            for (k, c0) in enumerate(qps)
                c1 = c0 + 1

                lon, lat = cell_centers(grid, c1)
                gc = golden["cell_centers"][k]
                @test close_rel(lon, gc["lon"], REL_TOL)
                @test close_rel(lat, gc["lat"], REL_TOL)

                got_nb = _neighbors_0idx(grid, c1)
                expected_nb = [Int(x) for x in golden["neighbors"][k]]
                @test got_nb == expected_nb

                for m in METRIC_NAMES
                    got = metric_eval(grid, m, c1)
                    expected = Float64(golden["metric_$(m)"][k])
                    @test close_rel(got, expected, REL_TOL)
                end
            end
        end
    end
end
