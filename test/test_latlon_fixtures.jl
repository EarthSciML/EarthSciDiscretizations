@testsnippet LatLonFixturesSetup begin
    using Test
    using JSON

    const FIXTURE_DIR = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "latlon",
    )
    const SCHEMA_PATH = joinpath(
        dirname(FIXTURE_DIR), "lat_lon.schema.json",
    )

    fixture_files() = sort([
        joinpath(FIXTURE_DIR, f)
            for f in readdir(FIXTURE_DIR)
            if endswith(f, ".esm")
    ])

    load_fixture(path) = JSON.parsefile(path)

    load_schema() = JSON.parsefile(SCHEMA_PATH)
end

@testitem "latlon fixtures: canonical set is present" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    files = fixture_files()
    names = [basename(f) for f in files]
    # Bead dsc-brl scope: 2-3 canonical resolutions (1°, 0.25°, 0.1°).
    @test "lat_lon_1deg.esm" in names
    @test "lat_lon_0p25deg.esm" in names
    @test "lat_lon_0p1deg.esm" in names
end

@testitem "latlon fixtures: declarative envelope matches GRIDS_API §7" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    for path in fixture_files()
        doc = load_fixture(path)
        # §7 common minimum + grid-level declarative fields.
        for key in ("family", "version", "dtype", "topology",
                "variant", "generator", "params")
            @test haskey(doc, key)
        end
        @test doc["family"] == "lat_lon"
        @test doc["version"] == "1.0.0"
        @test doc["topology"] == "rectilinear"
        @test doc["dtype"] in ("float64", "float32")
        @test doc["variant"] in ("regular", "reduced_gaussian")
        # Mayor 2026-04-20 correction: no inline geometry blobs in the fixture.
        for forbidden in ("cells", "edges", "vertices", "coordinates")
            @test !haskey(doc, forbidden)
        end
    end
end

@testitem "latlon fixtures: params conform to lat_lon.schema.json" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    schema = load_schema()
    @test schema["family"] == "lat_lon"
    known_opts = Set(keys(schema["options"]))

    for path in fixture_files()
        doc = load_fixture(path)
        params = doc["params"]
        @test params isa AbstractDict

        # Every params key must be declared in the family schema.
        for k in keys(params)
            @test k in known_opts
        end

        # Regular variant must carry nlon + nlat (schema §required).
        if doc["variant"] == "regular"
            @test haskey(params, "nlon")
            @test haskey(params, "nlat")
            @test params["nlon"] isa Integer && params["nlon"] >= 1
            @test params["nlat"] isa Integer && params["nlat"] >= 1
        end

        # Spherical family invariants: positive radius, non-negative halo,
        # only the implemented pole policy.
        R = get(params, "R", 6.371e6)
        @test R > 0
        @test get(params, "ghosts", 0) >= 0
        @test get(params, "pole_policy", "none") == "none"
    end
end

@testitem "lat_lon_1deg.esm pins expected canonical parameters" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    doc = load_fixture(joinpath(FIXTURE_DIR, "lat_lon_1deg.esm"))
    @test doc["generator"] == "lat_lon_regular"
    @test doc["params"]["nlon"] == 360
    @test doc["params"]["nlat"] == 180
    @test doc["params"]["R"] ≈ 6.371e6
    @test doc["params"]["ghosts"] == 0
    @test doc["params"]["lon_start"] ≈ -pi
end

@testitem "lat_lon_0p25deg.esm pins expected canonical parameters" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    doc = load_fixture(joinpath(FIXTURE_DIR, "lat_lon_0p25deg.esm"))
    @test doc["params"]["nlon"] == 1440
    @test doc["params"]["nlat"] == 720
    # 0.25° grid has exactly 4× the rows and columns of the 1° grid.
    @test doc["params"]["nlon"] == 4 * 360
    @test doc["params"]["nlat"] == 4 * 180
end

@testitem "lat_lon_0p1deg.esm pins expected canonical parameters" setup = [LatLonFixturesSetup] tags = [:grid, :lat_lon, :fixtures] begin
    doc = load_fixture(joinpath(FIXTURE_DIR, "lat_lon_0p1deg.esm"))
    @test doc["params"]["nlon"] == 3600
    @test doc["params"]["nlat"] == 1800
end
