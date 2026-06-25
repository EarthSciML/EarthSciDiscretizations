# Reprojection cross-binding conformance — REFERENCE BINDING (Julia) consumer.
#
# esd-47z.6 (RFC pure-io-data-loaders §8, "ESD DRAIN"). The reprojection/ family
# ships per-model Julia evaluator tests co-located with each .esm
# (reprojection/{longlat,lambert_conformal}_conformance_test.jl, beads esd-47z.1/.2).
# This file is the SHARED cross-binding layer: it consumes the binding-neutral
# fixtures + golden under tests/conformance/reprojection/<model>/ — the same
# fixtures.json/golden.json a Python/Rust/TypeScript consumer reads — and asserts
# the declarative documents are green THROUGH THE ESS ENGINE at the reference
# binding.
#
# The golden is an INDEPENDENT closed-form ORACLE (regenerate_golden.jl): longlat
# is the geographic identity (x = lon - lon_0); lambert_conformal is the spherical
# LCC Snyder closed form, cross-checked against proj4 reference points. Here the
# ESS engine evaluates the .esm via the zero-IC surfacing idiom (du = f!(u0) at the
# zero IC IS the surfaced map value) and must reproduce that oracle to the
# fixture's declared tolerance, plus the inverse-of-forward = identity round-trip.
# Per CONFORMANCE_SPEC §5.8 the non-Julia bindings consume the same golden; the
# cross-binding disposition is documented in tests/conformance/reprojection/README.md.

@testsnippet ReprojectionConformanceSetup begin
    using Test
    import EarthSciSerialization as ESS
    const JSON3 = ESS.JSON3
    using JSON   # the harness fixtures/golden are plain JSON

    const REPRO_HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "reprojection")
    const REPO_ROOT = joinpath(@__DIR__, "..")

    close_abs(a, b, atol) = abs(a - b) <= atol
    close_rel(a, b, rtol, atol) = abs(a - b) <= atol + rtol * abs(b)

    # Drive a reprojection .esm through build_evaluator at one coordinate, reading
    # the forward map (forward_x/forward_y) and the inverse-of-forward round-trip
    # (roundtrip_lon/roundtrip_lat) off the zero-IC surfacing states.
    function repro_eval(esm_relpath, model_name, overrides)
        raw = JSON3.read(read(joinpath(REPO_ROOT, esm_relpath), String))
        ics = Dict(
            "forward_x" => 0.0, "forward_y" => 0.0,
            "roundtrip_lon" => 0.0, "roundtrip_lat" => 0.0,
        )
        f!, u0, p, _, vmap = ESS.build_evaluator(
            raw; model_name = model_name,
            initial_conditions = ics, parameter_overrides = overrides,
        )
        du = similar(u0)
        f!(du, u0, p, 0.0)
        return (
            du[vmap["forward_x"]], du[vmap["forward_y"]],
            du[vmap["roundtrip_lon"]], du[vmap["roundtrip_lat"]],
        )
    end
end

@testitem "reprojection longlat cross-binding fixture (esd-47z.6)" setup = [ReprojectionConformanceSetup] tags = [:reprojection, :longlat, :conformance] begin
    dir = joinpath(REPRO_HARNESS, "longlat")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; model = fx["esm_model_name"]
    fwd_atol = Float64(fx["tolerance"]["forward_absolute"])
    rt_atol = Float64(fx["tolerance"]["roundtrip_absolute"])
    coords = fx["coords"]

    @test isfile(joinpath(REPO_ROOT, esm))
    @test length(gold["cases"]) == length(fx["cases"])

    for (ci, case) in enumerate(fx["cases"])
        gcase = gold["cases"][ci]
        @test gcase["name"] == case["name"]
        lon_0 = Float64(case["lon_0"])
        @testset "longlat case=$(case["name"])" begin
            for (k, c) in enumerate(coords)
                lon = Float64(c[1]); lat = Float64(c[2])
                fx_x, fx_y, rt_lon, rt_lat = repro_eval(
                    esm, model, Dict("lon" => lon, "lat" => lat, "lon_0" => lon_0),
                )
                gf = gcase["forward"][k]; gr = gcase["roundtrip"][k]
                # ENGINE reproduces the closed-form oracle (forward map)
                @test close_abs(fx_x, Float64(gf[1]), fwd_atol)
                @test close_abs(fx_y, Float64(gf[2]), fwd_atol)
                # inverse-of-forward = identity (round-trip), against golden == input
                @test close_abs(rt_lon, Float64(gr[1]), rt_atol)
                @test close_abs(rt_lat, Float64(gr[2]), rt_atol)
                @test close_abs(rt_lon, lon, rt_atol)
                @test close_abs(rt_lat, lat, rt_atol)
            end
        end
    end
end

@testitem "reprojection lambert_conformal cross-binding fixture (esd-47z.6)" setup = [ReprojectionConformanceSetup] tags = [:reprojection, :lambert_conformal, :conformance] begin
    dir = joinpath(REPRO_HARNESS, "lambert_conformal")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; model = fx["esm_model_name"]
    fwd_rtol = Float64(fx["tolerance"]["forward_relative"])
    fwd_atol = Float64(fx["tolerance"]["forward_absolute"])
    rt_atol = Float64(fx["tolerance"]["roundtrip_absolute"])
    coords = fx["coords"]

    @test isfile(joinpath(REPO_ROOT, esm))
    @test length(gold["param_sets"]) == length(fx["param_sets"])

    for (si, ps) in enumerate(fx["param_sets"])
        gps = gold["param_sets"][si]
        @test gps["name"] == ps["name"]
        p = ps["params"]
        @testset "lcc param_set=$(ps["name"])" begin
            for (k, c) in enumerate(coords)
                lon = Float64(c[1]); lat = Float64(c[2])
                overrides = Dict(
                    "lon" => lon, "lat" => lat,
                    "lat_1" => Float64(p["lat_1"]), "lat_2" => Float64(p["lat_2"]),
                    "lat_0" => Float64(p["lat_0"]), "lon_0" => Float64(p["lon_0"]),
                    "R" => Float64(p["R"]),
                )
                fx_x, fx_y, rt_lon, rt_lat = repro_eval(esm, model, overrides)
                gf = gps["forward"][k]
                # ENGINE reproduces the independent Snyder/proj4 oracle (forward map)
                @test close_rel(fx_x, Float64(gf[1]), fwd_rtol, fwd_atol)
                @test close_rel(fx_y, Float64(gf[2]), fwd_rtol, fwd_atol)
                # inverse-of-forward = identity (round-trip)
                @test close_abs(rt_lon, lon, rt_atol)
                @test close_abs(rt_lat, lat, rt_atol)
            end
            # false-origin invariant: the central-meridian point (k=1) has x == 0
            x0, _, _, _ = repro_eval(
                esm, model, Dict(
                    "lon" => Float64(coords[1][1]), "lat" => Float64(coords[1][2]),
                    "lat_1" => Float64(p["lat_1"]), "lat_2" => Float64(p["lat_2"]),
                    "lat_0" => Float64(p["lat_0"]), "lon_0" => Float64(p["lon_0"]),
                    "R" => Float64(p["R"]),
                ),
            )
            @test close_abs(x0, 0.0, 1.0e-6)
        end
    end
end
