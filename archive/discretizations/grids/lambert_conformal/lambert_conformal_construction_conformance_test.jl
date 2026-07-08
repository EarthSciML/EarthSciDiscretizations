# Projected-grid construction LAMBERT CONFORMAL — Julia evaluator conformance.
#
# ESD owns this declarative projected-grid construction (the grids/lambert_conformal/
# family: the wrf_lcc.esm / nei2016_lcc.esm Grid descriptors + rules/
# lambert_conformal_construction.esm); the EarthSciSerialization engine (an ESD
# dependency) EVALUATES it. This is the GDD Grid `crs` descriptor's first
# consumer (RFC pure-io-data-loaders §5.3, bead esd-47z.5): a PROJECTED native
# grid is topologically `cartesian` (a uniform lattice in projected metres)
# carrying a `crs` descriptor (projection=lambert_conformal + the WRF/NEI2016
# {lat_1, lat_2, lat_0, lon_0, R} parameters, ess-v9a.2); CONSTRUCTING its
# geographic geometry means recovering the (lon, lat) of every cell CORNER by
# applying the spherical-LCC INVERSE (rules/lambert_conformal_construction.esm,
# the inverse half of reprojection/lambert_conformal.esm, esd-47z.2) to the
# regular projected (x, y) corner lattice — closed-form, NO iteration:
#
#   ρ = sign(n)·sqrt(x² + (ρ₀−y)²)     θ = atan2(x, ρ₀−y)
#   lon = lon_0 + (θ/n)·180/π          lat = (2·atan((R·F/ρ)^(1/n)) − π/2)·180/π
#
# with the cone constants (n, F, R·F, ρ₀) fixed by the crs parameters. Like the
# rest of the reprojection/ family the construction is a SCALAR point-wise kernel
# (proj_x/proj_y scalar `parameter`s, geo_lon/geo_lat scalar `observed`s),
# surfaced through build_evaluator via the zero-IC state idiom: a state whose
# D-equation is d(state)/dt = <observed> from a zero IC, so du = f!(u0) IS the
# recovered coordinate. The harness here sweeps the kernel over the (Nx+1)×(Ny+1)
# corner lattice (a binding applies the point map over the corner array at the
# runtime layer).
#
# ACCEPTANCE (esd-47z.5): a projected GDD Grid round-trips and matches the
# WRF/NEI native geometry —
#   • Both projected Grid descriptors VALIDATE and load through the engine with
#     the `crs` descriptor preserved (the GDD Grid crs round-trip).
#   • The constructed cell-corner (lon, lat) match an INDEPENDENT proj4/Snyder
#     reference (the published spherical-LCC inverse, identical to PROJ
#     `+proj=lcc +R=…` for the sphere, computed in NumPy OUTSIDE the .esm) for
#     BOTH the WRF and NEI2016 grids — two real native grids, one construction.
#   • forward∘inverse = identity: forward-projecting the constructed corners
#     through the EXISTING reprojection/lambert_conformal.esm FORWARD map
#     restores the projected (x, y) lattice (the two reprojection-family members
#     compose to the identity — the grid round-trip).
#   • The crs parameters are LOAD-BEARING: WRF and NEI2016 (which share lon_0 =
#     −97 but differ in standard parallels, origin latitude, radius) recover
#     DIFFERENT geographic geometry from the IDENTICAL projected lattice off the
#     shared central meridian — the declarative-or-fail proof the construction
#     binds a genuinely parameterized grid.
# Per CONFORMANCE_SPEC §5.8 the other bindings schema-/structurally-validate the
# documents while Julia evaluates the construction here.

@testsnippet LCCGridConstructionSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach it
    # through the dependency rather than adding a redundant top-level dep — and so
    # the raw documents are the exact type build_evaluator expects.
    const JSON3 = ESS.JSON3

    const _CONSTRUCTION_FIXTURE = joinpath(@__DIR__, "rules", "lambert_conformal_construction.esm")
    const _WRF_GRID_FIXTURE = joinpath(@__DIR__, "wrf_lcc.esm")
    const _NEI2016_GRID_FIXTURE = joinpath(@__DIR__, "nei2016_lcc.esm")
    # The construction is the INVERSE half of the reprojection family; the
    # round-trip is closed by composing with the EXISTING forward map (top-level
    # reprojection/, three directories up from this grids/lambert_conformal/ dir).
    const _REPROJECTION_FIXTURE = joinpath(@__DIR__, "..", "..", "..", "reprojection", "lambert_conformal.esm")

    # INDEPENDENT inverse reference: (x, y, lon, lat) tuples computed from the
    # published spherical-LCC closed form (Snyder 1987 eqs. 15-5..15-9 = PROJ
    # `+proj=lcc +R=…` for the sphere) in NumPy, OUTSIDE the .esm AST, for the
    # representative CONUS lattice the two Grid descriptors pin (x0=−2e6, dx=1e6,
    # Nx=4; y0=−1.5e6, dy=1e6, Ny=3 ⇒ 5×4 corners). y outer, x inner (the corner
    # row-major order the harness builds from the fixture lattice parameters).
    const _GOLD_WRF = [
        (-2000000.0, -1500000.0, -116.10017801101773, 23.320805943967475),
        (-1000000.0, -1500000.0, -106.68789739792601, 24.883956953152033),
        (0.0, -1500000.0, -97.0, 25.415772540363083),
        (1000000.0, -1500000.0, -87.31210260207399, 24.883956953152033),
        (2000000.0, -1500000.0, -77.89982198898227, 23.320805943967475),
        (-2000000.0, -500000.0, -118.62443079701725, 31.92388539661814),
        (-1000000.0, -500000.0, -108.01300959338514, 33.76904204570805),
        (0.0, -500000.0, -97.0, 34.39829502737979),
        (1000000.0, -500000.0, -85.98699040661486, 33.76904204570805),
        (2000000.0, -500000.0, -75.37556920298275, 31.92388539661814),
        (-2000000.0, 500000.0, -121.89275561686046, 40.7292258345734),
        (-1000000.0, 500000.0, -109.75450863032499, 42.90012555470361),
        (0.0, 500000.0, -97.0, 43.64290686620447),
        (1000000.0, 500000.0, -84.24549136967501, 42.90012555470361),
        (2000000.0, 500000.0, -72.10724438313954, 40.7292258345734),
        (-2000000.0, 1500000.0, -126.27323349556194, 49.51039356753986),
        (-1000000.0, 1500000.0, -112.14243848502849, 52.06011386781491),
        (0.0, 1500000.0, -97.0, 52.93698916474947),
        (1000000.0, 1500000.0, -81.85756151497151, 52.06011386781491),
        (2000000.0, 1500000.0, -67.72676650443806, 49.51039356753986),
    ]
    const _GOLD_NEI2016 = [
        (-2000000.0, -1500000.0, -116.45427464977527, 24.64079229579825),
        (-1000000.0, -1500000.0, -106.83986695725679, 26.054764288944632),
        (0.0, -1500000.0, -97.0, 26.53329541521007),
        (1000000.0, -1500000.0, -87.16013304274321, 26.054764288944632),
        (2000000.0, -1500000.0, -77.54572535022473, 24.64079229579825),
        (-2000000.0, -500000.0, -118.7420023933866, 33.31438936166264),
        (-1000000.0, -500000.0, -108.0288124189089, 34.93283663909372),
        (0.0, -500000.0, -97.0, 35.4809903446274),
        (1000000.0, -500000.0, -85.9711875810911, 34.93283663909372),
        (2000000.0, -500000.0, -75.2579976066134, 33.31438936166264),
        (-2000000.0, 500000.0, -121.62456988567027, 42.05469408474961),
        (-1000000.0, 500000.0, -109.54249025768968, 43.890674083563056),
        (0.0, 500000.0, -97.0, 44.513136377440205),
        (1000000.0, 500000.0, -84.45750974231032, 43.890674083563056),
        (2000000.0, 500000.0, -72.37543011432973, 42.05469408474961),
        (-2000000.0, 1500000.0, -125.35974206311438, 50.63642210863358),
        (-1000000.0, 1500000.0, -111.53365889188319, 52.702256517538785),
        (0.0, 1500000.0, -97.0, 53.40378281423271),
        (1000000.0, 1500000.0, -82.46634110811681, 52.702256517538785),
        (2000000.0, 1500000.0, -68.64025793688562, 50.63642210863358),
    ]

    # Read a Grid descriptor's crs parameters from the fixture document (so the
    # fixtures are genuinely consumed — the kernel is driven by THEIR crs, and a
    # drift between fixture and reference would fail the native-geometry test).
    function _crs_params(grid_fixture::String, gridname::String)
        raw = JSON3.read(read(grid_fixture, String))
        crs = raw["grids"][gridname]["crs"]
        pr = crs["parameters"]
        return (
            lat_1 = Float64(pr["lat_1"]), lat_2 = Float64(pr["lat_2"]),
            lat_0 = Float64(pr["lat_0"]), lon_0 = Float64(pr["lon_0"]),
            R = Float64(crs["R"]),
        )
    end

    # Build the regular projected (x, y) corner lattice from the Grid descriptor's
    # own lattice parameters (x0, dx, Nx, y0, dy, Ny). (Nx+1)×(Ny+1) corner nodes
    # in y-outer/x-inner order, matching the reference tuple order.
    function _lattice_corners(grid_fixture::String, gridname::String)
        raw = JSON3.read(read(grid_fixture, String))
        pr = raw["grids"][gridname]["parameters"]
        Nx = Int(pr["Nx"]["default"]); Ny = Int(pr["Ny"]["default"])
        x0 = Float64(pr["x0"]["default"]); dx = Float64(pr["dx"]["default"])
        y0 = Float64(pr["y0"]["default"]); dy = Float64(pr["dy"]["default"])
        xs = [x0 + (i - 1) * dx for i in 1:(Nx + 1)]
        ys = [y0 + (j - 1) * dy for j in 1:(Ny + 1)]
        return [(x, y) for y in ys for x in xs]
    end

    # Apply the construction (spherical-LCC INVERSE) to a single projected corner
    # (proj_x, proj_y) under crs parameters p, returning the recovered
    # (corner_lon, corner_lat) in degrees. Each surfacing state is a constant-RHS
    # D-equation from a zero IC, so du = f!(u0) IS the recovered coordinate.
    function _construct_corner(x::Float64, y::Float64, p)
        raw = JSON3.read(read(_CONSTRUCTION_FIXTURE, String))
        ics = Dict("corner_lon" => 0.0, "corner_lat" => 0.0)
        f!, u0, pp, _, vmap = ESS.build_evaluator(
            raw; model_name = "LambertConformalGridConstruction",
            initial_conditions = ics,
            parameter_overrides = Dict(
                "proj_x" => x, "proj_y" => y,
                "lat_1" => p.lat_1, "lat_2" => p.lat_2, "lat_0" => p.lat_0,
                "lon_0" => p.lon_0, "R" => p.R,
            ),
        )
        du = similar(u0)
        f!(du, u0, pp, 0.0)
        return (du[vmap["corner_lon"]], du[vmap["corner_lat"]])
    end

    # Forward-project a geographic (lon, lat) through the EXISTING reprojection
    # rule's FORWARD map under crs parameters p, returning (forward_x, forward_y)
    # in metres — the other half of the round-trip.
    function _reproject_forward(lon::Float64, lat::Float64, p)
        raw = JSON3.read(read(_REPROJECTION_FIXTURE, String))
        ics = Dict(
            "forward_x" => 0.0, "forward_y" => 0.0,
            "roundtrip_lon" => 0.0, "roundtrip_lat" => 0.0,
        )
        f!, u0, pp, _, vmap = ESS.build_evaluator(
            raw; model_name = "LambertConformal",
            initial_conditions = ics,
            parameter_overrides = Dict(
                "lon" => lon, "lat" => lat,
                "lat_1" => p.lat_1, "lat_2" => p.lat_2, "lat_0" => p.lat_0,
                "lon_0" => p.lon_0, "R" => p.R,
            ),
        )
        du = similar(u0)
        f!(du, u0, pp, 0.0)
        return (du[vmap["forward_x"]], du[vmap["forward_y"]])
    end
end

@testitem "projected grid lambert_conformal construction (esd-47z.5)" setup = [LCCGridConstructionSetup] tags = [:grids, :lambert_conformal, :crs, :reprojection, :conformance] begin

    # Both projected Grid descriptors VALIDATE and load through the engine with
    # the GDD Grid `crs` descriptor (ess-v9a.2). This is the crs round-trip at the
    # document level: a projected native grid (WRF/NEI LCC params) parses, and the
    # crs survives as a lambert_conformal descriptor with the expected projection.
    @testset "projected Grid descriptors load with crs (schema + structural)" begin
        for (fixture, gridname) in (
                (_WRF_GRID_FIXTURE, "wrf_lcc"), (_NEI2016_GRID_FIXTURE, "nei2016_lcc"),
            )
            @test isfile(fixture)
            @test (ESS.load(fixture); true)
            raw = JSON3.read(read(fixture, String))
            crs = raw["grids"][gridname]["crs"]
            @test String(crs["projection"]) == "lambert_conformal"
            @test String(crs["datum"]) == "sphere"
            @test raw["grids"][gridname]["family"] == "cartesian"  # crs is orthogonal to family
        end
        @test isfile(_CONSTRUCTION_FIXTURE)
        @test (ESS.load(_CONSTRUCTION_FIXTURE); true)
    end

    # The harness builds the projected lattice from each fixture's OWN lattice
    # parameters; confirm the (Nx+1)×(Ny+1) corner set is exactly the reference
    # lattice (so the fixture's lattice spec is what's under test, and the corner
    # order matches the reference tuples).
    @testset "fixture lattice parameters reproduce the reference corner lattice" begin
        for (fixture, gridname, gold) in (
                (_WRF_GRID_FIXTURE, "wrf_lcc", _GOLD_WRF),
                (_NEI2016_GRID_FIXTURE, "nei2016_lcc", _GOLD_NEI2016),
            )
            corners = _lattice_corners(fixture, gridname)
            @test length(corners) == length(gold)
            for (k, (x, y)) in enumerate(corners)
                @test x == gold[k][1]
                @test y == gold[k][2]
            end
        end
    end

    # ACCEPTANCE — the construction matches the WRF/NEI native geometry. The .esm
    # AST and the NumPy reference both evaluate the SAME closed-form spherical-LCC
    # inverse, differing only by cross-libm rounding over the ~12-op chain, so they
    # agree to ≪ µdeg; a transcription error in the AST would be whole degrees off,
    # far outside the tolerance. Two real native grids exercise the ONE construction.
    @testset "construction matches proj4/Snyder native geometry (WRF + NEI2016)" begin
        for (fixture, gridname, gold) in (
                (_WRF_GRID_FIXTURE, "wrf_lcc", _GOLD_WRF),
                (_NEI2016_GRID_FIXTURE, "nei2016_lcc", _GOLD_NEI2016),
            )
            p = _crs_params(fixture, gridname)
            for (x, y, glon, glat) in gold
                clon, clat = _construct_corner(x, y, p)
                @test isapprox(clon, glon; rtol = 1.0e-9, atol = 1.0e-8)
                @test isapprox(clat, glat; rtol = 1.0e-9, atol = 1.0e-8)
            end
        end
    end

    # ACCEPTANCE — forward∘inverse = identity (the grid round-trip). Forward-
    # projecting each constructed corner through the EXISTING reprojection forward
    # map restores the projected (x, y) lattice node: the construction (inverse)
    # and reprojection (forward) share the cone constants by construction, so they
    # compose to the identity. On the shared central meridian (x = 0) the forward
    # easting is exactly 0 (θ = 0).
    @testset "forward∘inverse restores the projected lattice (round-trip)" begin
        for (fixture, gridname, gold) in (
                (_WRF_GRID_FIXTURE, "wrf_lcc", _GOLD_WRF),
                (_NEI2016_GRID_FIXTURE, "nei2016_lcc", _GOLD_NEI2016),
            )
            p = _crs_params(fixture, gridname)
            for (x, y, _glon, _glat) in gold
                clon, clat = _construct_corner(x, y, p)
                fx, fy = _reproject_forward(clon, clat, p)
                @test isapprox(fx, x; rtol = 1.0e-7, atol = 1.0e-3)
                @test isapprox(fy, y; rtol = 1.0e-7, atol = 1.0e-3)
            end
        end
    end

    # The crs parameters are LOAD-BEARING (the declarative-or-fail proof the
    # construction binds a PARAMETERIZED grid, not a hard-coded one): WRF and
    # NEI2016 share lon_0 = −97 but differ in standard parallels (cone constant
    # n ≈ 0.7156 vs 0.6305), origin latitude, and radius, so the IDENTICAL
    # projected lattice recovers DIFFERENT geographic geometry off the shared
    # central meridian. On the central meridian (x = 0) both recover lon = −97
    # exactly (θ = 0), so that column is excluded from the inequality.
    @testset "crs parameters are load-bearing (WRF ≠ NEI2016)" begin
        pW = _crs_params(_WRF_GRID_FIXTURE, "wrf_lcc")
        pN = _crs_params(_NEI2016_GRID_FIXTURE, "nei2016_lcc")
        for (x, y, _glon, _glat) in _GOLD_WRF
            lonW, latW = _construct_corner(x, y, pW)
            lonN, latN = _construct_corner(x, y, pN)
            if x == 0.0
                @test isapprox(lonW, -97.0; atol = 1.0e-9)  # shared central meridian
                @test isapprox(lonN, -97.0; atol = 1.0e-9)
            else
                @test lonW != lonN
            end
            @test latW != latN  # latitude differs at every corner (different n, lat_0, R)
        end
    end
end
