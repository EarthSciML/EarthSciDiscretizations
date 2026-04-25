using Test
using TestItems

# Tests for the three-layer CI walker that validates rule files under
# discretizations/. The walker discovers rules via load_rules and runs three
# layers per rule: (A) canonical-form byte-diff, (B) MMS convergence (driven
# by EarthSciSerialization.verify_mms_convergence), (C) integration
# benchmarks. Layer A skips until canonical-form fixtures are authored;
# Layer C is gated on ESD_RUN_INTEGRATION=1 and skipped by default.

@testitem "walker: discovers seeded rules; centered_2nd_uniform layer B passes via ESS evaluator" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    junit = joinpath(@__DIR__, "junit-esd.xml")
    isfile(junit) && rm(junit)

    # Pin Layer C to the default skip path so the broad assertion below holds
    # regardless of whether the caller exported ESD_RUN_INTEGRATION=1.
    prior_run_integration = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    delete!(ENV, "ESD_RUN_INTEGRATION")

    results = try
        WalkESDTests.walk_esd_tests(; catalog = catalog, junit_path = junit)
    finally
        if prior_run_integration === nothing
            delete!(ENV, "ESD_RUN_INTEGRATION")
        else
            ENV["ESD_RUN_INTEGRATION"] = prior_run_integration
        end
    end

    names = Set(r.name for r in results)
    # Superset: the catalog has grown with grid schemas and other families;
    # only require that the canonical FD rules remain walker-discovered.
    # `centered_2nd_uniform_vertical` (dsc-mzu) is the first non-cartesian
    # rule and a design probe for the per-family selector schema.
    for seeded in (
        "centered_2nd_uniform",
        "centered_2nd_uniform_vertical",
        "centered_2nd_uniform_latlon",
        "nn_diffusion_mpas",
        "periodic_bc",
        "upwind_1st",
    )
        @test seeded in names
    end
    # finite_volume/ppm_reconstruction (CW84 §1) joins the walker once landed.
    @test "ppm_reconstruction" in names
    ppm = first(filter(r -> r.name == "ppm_reconstruction", results))
    @test ppm.family == :finite_volume
    # finite_volume/weno5_advection (Jiang-Shu 1996) joins the walker too.
    @test "weno5_advection" in names
    weno = first(filter(r -> r.name == "weno5_advection", results))
    @test weno.family == :finite_volume

    # Layer A always skips (ESS canonical-form rule engine not yet wired).
    # Layer C always skips unless ESD_RUN_INTEGRATION=1.
    # Layer B passes for rules with a runnable convergence fixture
    # (centered_2nd_uniform, centered_2nd_uniform_vertical, upwind_1st —
    # linear stencils evaluated via the ESS AST evaluator). Rules whose
    # convergence fixture declares applicable:false (limiters, reconstruction
    # rules pending ESS harness extension, periodic_bc) skip with a
    # fixture-declared reason. Rules with no convergence fixture at all also
    # skip.
    pass_layer_b = Set([("finite_difference", "centered_2nd_uniform"),
                        ("finite_difference", "centered_2nd_uniform_vertical"),
                        ("finite_difference", "upwind_1st")])
    not_applicable_layer_b = Set([("finite_difference", "periodic_bc"),
                                   ("finite_difference", "centered_2nd_uniform_latlon"),
                                   ("finite_difference", "nn_diffusion_mpas"),
                                   ("finite_volume", "ppm_reconstruction"),
                                   ("finite_volume", "weno5_advection"),
                                   ("finite_volume", "flux_limiter_minmod"),
                                   ("finite_volume", "flux_limiter_superbee")])
    for r in results
        @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
        @test !isempty(r.layer_a.reason)
        @test r.layer_c.outcome == WalkESDTests.LAYER_SKIP
        @test !isempty(r.layer_c.reason)
        key = (String(r.family), r.name)
        if key in pass_layer_b
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif key in not_applicable_layer_b
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        else
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test !isempty(r.layer_b.reason)
        end
    end

    @test isfile(junit)
    xml = read(junit, String)
    @test occursin("<testsuites>", xml)
    @test occursin("<testsuite name=\"ESD Walker\"", xml)
    # Parametrize against actual catalog size: 3 layers (A/B/C) per rule.
    total = length(results) * 3
    # Three layer-B cases pass (centered_2nd_uniform,
    # centered_2nd_uniform_vertical, upwind_1st); the rest skip.
    passed = sum(1 for r in results
                 if (String(r.family), r.name) in pass_layer_b; init = 0)
    skipped = total - passed
    @test passed == 3
    @test occursin("tests=\"$total\"", xml)
    @test occursin("failures=\"0\"", xml)
    @test occursin("skipped=\"$skipped\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_vertical\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_latlon\"", xml)
    @test occursin("classname=\"finite_volume.ppm_reconstruction\"", xml)
    @test occursin("classname=\"finite_volume.weno5_advection\"", xml)
    @test occursin("name=\"layer_A\"", xml)
    @test occursin("<skipped message=", xml)
end

@testitem "walker: empty catalog yields zero results and a green report" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests

    mktempdir() do tmp
        # Empty directory: no families, no rules.
        results = WalkESDTests.walk_esd_tests(; catalog = tmp)
        @test isempty(results)
    end

    # Non-existent path is treated the same way (zero rules, no error).
    @test isempty(WalkESDTests.discover_rules(joinpath(tempdir(), "nope-$(rand(UInt32))")))
end

@testitem "walker: layer A flags missing fixture files as failure" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_difference")
        mkpath(family_dir)
        rule_json = joinpath(family_dir, "broken_rule.json")
        write(rule_json, "{}")
        # Create an empty canonical fixture directory: walker should detect
        # the directory is present but inputs are missing and surface a fail.
        canonical = joinpath(family_dir, "broken_rule", "fixtures", "canonical")
        mkpath(canonical)

        rule = RuleFile(:finite_difference, "broken_rule", rule_json)
        result = WalkESDTests.run_layer_a(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("missing input.esm or expected.esm", result.reason)
    end
end

@testitem "walker: layer C honors ESD_RUN_INTEGRATION env var" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_difference")
        mkpath(family_dir)
        rule_json = joinpath(family_dir, "any_rule.json")
        write(rule_json, "{}")
        rule = RuleFile(:finite_difference, "any_rule", rule_json)

        prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
        try
            delete!(ENV, "ESD_RUN_INTEGRATION")
            r = WalkESDTests.run_layer_c(rule)
            @test r.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("ESD_RUN_INTEGRATION=1", r.reason)

            ENV["ESD_RUN_INTEGRATION"] = "1"
            r2 = WalkESDTests.run_layer_c(rule)
            @test r2.outcome == WalkESDTests.LAYER_SKIP
            # No integration fixtures exist for this synthetic rule, so the
            # skip reason switches to "no integration fixtures declared".
            @test occursin("no integration fixtures", r2.reason)
        finally
            if prior === nothing
                delete!(ENV, "ESD_RUN_INTEGRATION")
            else
                ENV["ESD_RUN_INTEGRATION"] = prior
            end
        end
    end
end

@testitem "walker: layer C runs Cartesian Gaussian end-to-end via ESS pipeline (ESD_RUN_INTEGRATION=1)" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    upwind = first(filter(r -> r.name == "upwind_1st", rules))

    prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    try
        ENV["ESD_RUN_INTEGRATION"] = "1"
        result = WalkESDTests.run_layer_c(upwind)
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("gaussian_advection_cartesian_1d", result.reason)
        @test occursin("L∞", result.reason)
    finally
        if prior === nothing
            delete!(ENV, "ESD_RUN_INTEGRATION")
        else
            ENV["ESD_RUN_INTEGRATION"] = prior
        end
    end
end

@testitem "walker: layer C runs Williamson 1 cubed-sphere advection (ESD_RUN_INTEGRATION=1)" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    ppm = first(filter(r -> r.name == "ppm_reconstruction", rules))

    prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    try
        ENV["ESD_RUN_INTEGRATION"] = "1"
        result = WalkESDTests.run_layer_c(ppm)
        # Williamson 1 PASSes; Williamson 2 + DCMIP 1-1 are stubs that SKIP.
        # Aggregate: 1 pass + 2 skip → LAYER_PASS with "1/3 cases pass, 2 skipped".
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("williamson1_cosine_bell", result.reason)
        @test occursin("williamson2_geostrophic_steady", result.reason)
        @test occursin("dcmip_1_1_3d_advection", result.reason)
    finally
        if prior === nothing
            delete!(ENV, "ESD_RUN_INTEGRATION")
        else
            ENV["ESD_RUN_INTEGRATION"] = prior
        end
    end
end

@testitem "walker: junit XML escapes special characters" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests

    results = [
        WalkESDTests.RuleResult(
            :finite_difference,
            "rule_with_<>&\"'",
            "/tmp/x.json",
            WalkESDTests.LayerResult(WalkESDTests.LAYER_SKIP, "reason with <tag> & \"quote\""),
            WalkESDTests.LayerResult(WalkESDTests.LAYER_SKIP, ""),
            WalkESDTests.LayerResult(WalkESDTests.LAYER_SKIP, ""),
        ),
    ]
    mktempdir() do tmp
        path = joinpath(tmp, "junit.xml")
        WalkESDTests.write_junit(path, results)
        xml = read(path, String)
        @test occursin("&lt;tag&gt;", xml)
        @test occursin("&amp;", xml)
        @test occursin("&quot;quote&quot;", xml)
        @test !occursin("<tag>", xml)
    end
end
