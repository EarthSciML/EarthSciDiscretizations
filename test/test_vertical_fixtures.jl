@testsnippet VerticalFixturesSetup begin
    using Test
    using JSON

    const VERTICAL_FIXTURE_DIR = joinpath(
        @__DIR__, "..", "discretizations", "grids", "vertical",
    )

    const VERTICAL_SCHEMA_PATH = joinpath(
        @__DIR__, "..", "discretizations", "grids", "vertical.schema.json",
    )

    # Coordinate kinds the vertical family accepts. Must stay in sync with
    # `python/src/earthsci_toolkit/grids/vertical.py::_COORDINATES` and its
    # Julia / Rust / TS siblings.
    const VERTICAL_COORDINATES = Set(
        [
            "sigma", "eta", "z", "theta", "hybrid_sigma_theta", "z_star",
        ]
    )

    # Coordinates whose `levels` are sigma-like (strictly decreasing, 1 →
    # surface, 0 → top). The remaining coordinates have strictly-increasing
    # `levels` in their native units.
    const SIGMA_LIKE_COORDS = Set(["sigma", "eta", "hybrid_sigma_theta"])

    function list_esm_fixtures(dir::AbstractString)
        isdir(dir) || return String[]
        return sort(
            [
                joinpath(dir, f) for f in readdir(dir) if endswith(f, ".esm")
            ]
        )
    end
end

@testitem "vertical fixtures: schema file is valid" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    @test isfile(VERTICAL_SCHEMA_PATH)
    schema = JSON.parsefile(VERTICAL_SCHEMA_PATH)
    @test schema["family"] == "vertical"
    @test schema["version"] == "1.0.0"
    @test "coordinate" in schema["required"]
    @test haskey(schema["options"], "coordinate")
    @test haskey(schema["options"], "levels")
    @test haskey(schema["options"], "ak")
    @test haskey(schema["options"], "bk")
    @test haskey(schema["options"], "p0")
    @test haskey(schema["options"], "dtype")
    @test haskey(schema["options"], "ghosts")
    @test schema["tolerances"]["default_rel"] ≈ 1.0e-14
end

@testitem "vertical fixtures: directory contains .esm files" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    fixtures = list_esm_fixtures(VERTICAL_FIXTURE_DIR)
    @test !isempty(fixtures)
    # Must carry at least one fixture per bead acceptance (2-3 canonical
    # resolutions). The per-coordinate coverage below is asserted separately.
    @test length(fixtures) >= 2
end

@testitem "vertical fixtures: every .esm parses and has the common-minimum shape" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    fixtures = list_esm_fixtures(VERTICAL_FIXTURE_DIR)
    @test !isempty(fixtures)
    for path in fixtures
        @testset "$(basename(path))" begin
            doc = JSON.parsefile(path)
            # GRIDS_API.md §7 common-minimum fields
            @test doc["family"] == "vertical"
            @test doc["topology"] == "column"
            @test doc["ndim"] == 1
            @test doc["dtype"] in ("float64", "float32")
            @test doc["ghosts"] isa Integer
            @test doc["ghosts"] >= 0
            @test doc["n_cells"] isa Integer
            @test doc["n_vertices"] == doc["n_cells"] + 1
            @test doc["n_edges"] == doc["n_cells"]
            @test doc["schema_version"] == "1.0.0"

            # Declarative options block
            opts = doc["options"]
            @test opts["coordinate"] in VERTICAL_COORDINATES
            @test opts["nz"] == doc["n_cells"]
            levels = opts["levels"]
            @test levels isa AbstractVector
            @test length(levels) == opts["nz"] + 1

            # Provenance block
            prov = doc["provenance"]
            @test haskey(prov, "binding")
            @test prov["family"] == "vertical"
            @test prov["version"] == "1.0.0"
            @test prov["coordinate"] == opts["coordinate"]
            @test prov["dtype"] == doc["dtype"]
        end
    end
end

@testitem "vertical fixtures: levels respect the coordinate monotonicity contract" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    for path in list_esm_fixtures(VERTICAL_FIXTURE_DIR)
        @testset "$(basename(path))" begin
            doc = JSON.parsefile(path)
            coord = doc["options"]["coordinate"]
            levels = Float64.(doc["options"]["levels"])
            diffs = diff(levels)
            if coord in SIGMA_LIKE_COORDS
                @test all(d -> d < 0, diffs)  # strictly decreasing
                @test levels[1] ≈ 1.0
                @test abs(levels[end]) < 1.0e-12
            else
                @test all(d -> d > 0, diffs)  # strictly increasing
            end
        end
    end
end

@testitem "vertical fixtures: hybrid coefficients are consistent where present" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    for path in list_esm_fixtures(VERTICAL_FIXTURE_DIR)
        @testset "$(basename(path))" begin
            doc = JSON.parsefile(path)
            opts = doc["options"]
            coord = opts["coordinate"]
            nz = opts["nz"]
            if coord == "eta"
                @test haskey(opts, "ak")
                @test haskey(opts, "bk")
                @test haskey(opts, "p0")
                ak = Float64.(opts["ak"])
                bk = Float64.(opts["bk"])
                @test length(ak) == nz + 1
                @test length(bk) == nz + 1
                p0 = Float64(opts["p0"])
                @test p0 > 0
                @test isfinite(p0)
                # Synthesized sigma = ak/p0 + bk must be strictly decreasing,
                # with eta's own `levels` agreeing with it (the Python
                # binding stores the synthesized sigma in `levels`).
                sigma = ak ./ p0 .+ bk
                @test all(d -> d < 0, diff(sigma))
                @test isapprox(sigma, Float64.(opts["levels"]); atol = 1.0e-12)
            elseif coord == "hybrid_sigma_theta"
                # ak/bk are optional here. When present they must match levels.
                if haskey(opts, "ak")
                    @test length(opts["ak"]) == nz + 1
                end
                if haskey(opts, "bk")
                    @test length(opts["bk"]) == nz + 1
                end
            else
                @test !haskey(opts, "ak")
                @test !haskey(opts, "bk")
            end
        end
    end
end

@testitem "vertical fixtures: canonical resolutions cover the required coordinate kinds" setup = [VerticalFixturesSetup] tags = [:grid, :vertical, :fixtures] begin
    coords_seen = Set{String}()
    nz_by_coord = Dict{String, Vector{Int}}()
    for path in list_esm_fixtures(VERTICAL_FIXTURE_DIR)
        doc = JSON.parsefile(path)
        coord = doc["options"]["coordinate"]
        push!(coords_seen, coord)
        push!(get!(nz_by_coord, coord, Int[]), Int(doc["options"]["nz"]))
    end
    # Sigma is the coarsest coordinate; require the two-resolution bead
    # acceptance there. Other coordinate kinds contribute at least one
    # fixture each so the family's declarative schema is exercised end-to-end.
    @test "sigma" in coords_seen
    @test length(nz_by_coord["sigma"]) >= 2
    @test "z" in coords_seen
    @test "eta" in coords_seen
    @test "theta" in coords_seen
end
