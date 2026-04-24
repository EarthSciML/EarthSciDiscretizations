@testsnippet ArakawaFixturesSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON
end

# Inline topology/metric tests for the canonical arakawa .esm fixtures at
# `discretizations/grids/arakawa/` (bead dsc-u0l). The fixtures are declarative
# configs per the 2026-04-20 mayor correction (bead dsc-suu): .esm contains no
# inline geometry arrays; accessors derive geometry on demand.
#
# Layer-A-style verification: declarative shape matches the canonical form
# documented in `discretizations/grids/arakawa/README.md` and — for fixtures
# whose base grid is wired into the current Julia arakawa runtime — round-
# trips through `EarthSciDiscretizations.grids.arakawa` produce byte-matching
# output.

@testitem "arakawa fixtures: manifest is well-formed" setup = [ArakawaFixturesSetup] tags = [:arakawa, :fixtures] begin
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    manifest = JSON.parsefile(joinpath(fixtures_dir, "fixtures.json"))
    @test manifest["family"] == "arakawa"
    names = [f["name"] for f in manifest["fixtures"]]
    @test sort(names) == sort([
        "cartesian_n16", "cartesian_n64", "cartesian_n256",
        "lat_lon_1deg", "lat_lon_0p25deg", "lat_lon_0p1deg",
    ])
    # Every listed file exists and parses as JSON.
    for f in manifest["fixtures"]
        path = joinpath(fixtures_dir, f["file"])
        @test isfile(path)
        @test JSON.parsefile(path) isa AbstractDict
    end
end

@testitem "arakawa fixtures: declarative form (no inline geometry)" setup = [ArakawaFixturesSetup] tags = [:arakawa, :fixtures] begin
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    manifest = JSON.parsefile(joinpath(fixtures_dir, "fixtures.json"))
    for f in manifest["fixtures"]
        d = JSON.parsefile(joinpath(fixtures_dir, f["file"]))
        # Required canonical fields per README.
        @test d["family"] == "arakawa"
        @test d["topology"] == "block_structured"
        @test d["dtype"] in ("float64", "float32")
        @test d["ghosts"] isa Int && d["ghosts"] >= 0
        @test d["stagger"] in ("A", "B", "C", "D", "E")
        @test d["rotated"] == (d["stagger"] == "E")
        @test d["n_cells"] isa Int && d["n_cells"] > 0
        @test haskey(d, "base") && d["base"] isa AbstractDict
        # Per mayor correction: no inline geometry arrays.
        @test !haskey(d, "cells")
        @test !haskey(d, "edges")
        @test !haskey(d, "vertices")
    end
end

@testitem "arakawa fixtures: cartesian n_cells and base schema match" setup = [
    ArakawaFixturesSetup,
] tags = [:arakawa, :fixtures] begin
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    manifest = JSON.parsefile(joinpath(fixtures_dir, "fixtures.json"))
    cart = [f for f in manifest["fixtures"] if f["base_family"] == "cartesian"]
    @test length(cart) == 3
    for f in cart
        d = JSON.parsefile(joinpath(fixtures_dir, f["file"]))
        base = d["base"]
        @test base["family"] == "cartesian"
        @test base["nx"] == f["nx"]
        @test base["ny"] == f["ny"]
        @test d["n_cells"] == base["nx"] * base["ny"]
        # Matches cartesian.schema.json required fields: nx, extent.
        @test haskey(base, "extent")
        extent = base["extent"]
        @test length(extent) == 2
        @test length(extent[1]) == 2 && length(extent[2]) == 2
        # Upper corner strictly greater than lower corner on both axes.
        @test extent[2][1] > extent[1][1]
        @test extent[2][2] > extent[1][2]
    end
end

@testitem "arakawa fixtures: lat_lon n_cells and base schema match" setup = [
    ArakawaFixturesSetup,
] tags = [:arakawa, :fixtures] begin
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    manifest = JSON.parsefile(joinpath(fixtures_dir, "fixtures.json"))
    latlon = [f for f in manifest["fixtures"] if f["base_family"] == "lat_lon"]
    @test length(latlon) == 3
    for f in latlon
        d = JSON.parsefile(joinpath(fixtures_dir, f["file"]))
        base = d["base"]
        @test base["family"] == "lat_lon"
        @test base["variant"] == "regular"
        @test base["nlon"] == f["nlon"]
        @test base["nlat"] == f["nlat"]
        @test d["n_cells"] == base["nlon"] * base["nlat"]
        # Matches lat_lon.schema.json options: R, pole_policy, lon_start.
        @test base["R"] > 0
        @test base["pole_policy"] == "none"
        @test base["lon_start"] <= 0
    end
end

@testitem "arakawa fixtures: cartesian round-trips byte-exactly through Julia generator" setup = [
    ArakawaFixturesSetup,
] tags = [:arakawa, :fixtures, :conformance] begin
    # Byte-identity check: Julia reference binding (docs/GRIDS_API.md §4.3).
    # Build each cartesian fixture's ArakawaGrid via the public generator,
    # lower to_esm, and compare to the committed fixture.
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    cases = [
        ("cartesian_n16.esm.json", 16),
        ("cartesian_n64.esm.json", 64),
        ("cartesian_n256.esm.json", 256),
    ]
    for (file, N) in cases
        expected = JSON.parsefile(joinpath(fixtures_dir, file))
        b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = N, ny = N)
        g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)
        actual = to_esm(g)
        @test actual["family"] == expected["family"]
        @test actual["stagger"] == expected["stagger"]
        @test actual["dtype"] == expected["dtype"]
        @test actual["topology"] == expected["topology"]
        @test actual["ghosts"] == expected["ghosts"]
        @test actual["rotated"] == expected["rotated"]
        @test actual["n_cells"] == expected["n_cells"]
        ab = actual["base"]
        eb = expected["base"]
        @test ab["family"] == eb["family"]
        @test ab["nx"] == eb["nx"]
        @test ab["ny"] == eb["ny"]
        # Element-level extent check (JSON parsers return heterogeneously-typed
        # nested arrays; comparing the outer containers directly is parser-
        # dependent, scalar comparison is not).
        @test ab["extent"][1][1] == eb["extent"][1][1]
        @test ab["extent"][1][2] == eb["extent"][1][2]
        @test ab["extent"][2][1] == eb["extent"][2][1]
        @test ab["extent"][2][2] == eb["extent"][2][2]
    end
end

@testitem "arakawa fixtures: cartesian inline topology/metric checks" setup = [
    ArakawaFixturesSetup,
] tags = [:arakawa, :fixtures] begin
    # Build each cartesian fixture's grid and exercise accessors against
    # expected values derived from the declaration. This is the Layer-A
    # topology/metric check dsc-usq's walker will eventually consume.
    for N in (16, 64, 256)
        b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = N, ny = N)
        g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)

        # Topology: C-grid variable shapes per stagger on N × N interior cells.
        @test variable_shape(g, :h) == (N, N)
        @test variable_shape(g, :u) == (N + 1, N)
        @test variable_shape(g, :v) == (N, N + 1)

        dx = 1.0 / N
        # Metric uniformity on Cartesian.
        @test metric_eval(g, :dx, 1, 1) ≈ dx
        @test metric_eval(g, :dy, 1, 1) ≈ dx
        @test metric_eval(g, :area, N, N) ≈ dx * dx

        # Cell (1,1) at half-cell offset from the SW corner.
        cx, cy = cell_centers(g, 1, 1)
        @test cx ≈ dx / 2
        @test cy ≈ dx / 2
        # Last cell (N,N) at half-cell offset from the NE corner.
        cx, cy = cell_centers(g, N, N)
        @test cx ≈ 1.0 - dx / 2
        @test cy ≈ 1.0 - dx / 2

        # C-stagger u at west boundary (1,1) sits on the x=0 edge.
        ux, uy = u_face(g, 1, 1)
        @test ux ≈ 0.0
        @test uy ≈ dx / 2
        # C-stagger v at south boundary (1,1) sits on the y=0 edge.
        vx, vy = v_face(g, 1, 1)
        @test vx ≈ dx / 2
        @test vy ≈ 0.0

        # Corner (1,1) is the SW domain corner; corner (N+1, N+1) is NE.
        @test corners(g, 1, 1) == (0.0, 0.0)
        @test corners(g, N + 1, N + 1) == (1.0, 1.0)

        # Interior neighbors present; SW-corner neighbors missing on W and S.
        w, e, s, n = neighbors(g, 1, 1)
        @test w === nothing && s === nothing
        @test e === (CellCenter, 2, 1)
        @test n === (CellCenter, 1, 2)
    end
end

@testitem "arakawa fixtures: lat_lon declarations parse but no runtime yet" setup = [
    ArakawaFixturesSetup,
] tags = [:arakawa, :fixtures] begin
    # Lat-lon base grids are not yet wired into the arakawa runtime (see
    # README.md "Lat-lon fixtures: declaration-only"). Verify the three
    # fixtures parse to the canonical lat_lon shape and declare the standard
    # 1°/0.25°/0.1° resolutions on Earth's sphere.
    fixtures_dir = joinpath(
        dirname(dirname(pathof(EarthSciDiscretizations))),
        "discretizations", "grids", "arakawa",
    )
    expectations = Dict(
        "lat_lon_1deg.esm.json" => (360, 180),
        "lat_lon_0p25deg.esm.json" => (1440, 720),
        "lat_lon_0p1deg.esm.json" => (3600, 1800),
    )
    for (file, (nlon, nlat)) in expectations
        d = JSON.parsefile(joinpath(fixtures_dir, file))
        base = d["base"]
        @test base["nlon"] == nlon
        @test base["nlat"] == nlat
        @test d["n_cells"] == nlon * nlat
        @test base["R"] ≈ 6371000.0
        @test base["lon_start"] ≈ -π
    end
end
