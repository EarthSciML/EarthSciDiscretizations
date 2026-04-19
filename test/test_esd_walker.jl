using Test
using TestItems

# Tests for the three-layer CI walker that validates rule files under
# discretizations/. The walker discovers rules via load_rules and runs three
# layers per rule: (A) canonical-form byte-diff, (B) MMS convergence, (C)
# integration benchmarks. Layers A and B require dependencies that have not
# yet landed (ESS rule engine gt-b13f and the tree-walk evaluator), so for
# now they skip with descriptive reasons. Layer C is gated on
# ESD_RUN_INTEGRATION=1 and skipped by default.

@testitem "walker: discovers the three seeded rules and skips with reasons" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    junit = joinpath(@__DIR__, "junit-esd.xml")
    isfile(junit) && rm(junit)

    results = WalkESDTests.walk_esd_tests(; catalog = catalog, junit_path = junit)

    names = sort([r.name for r in results])
    @test names == ["centered_2nd_uniform", "periodic_bc", "upwind_1st"]

    # No fixtures authored yet, so every layer should skip cleanly. The
    # reason string must be non-empty so the JUnit consumer surfaces it.
    for r in results
        for layer in (r.layer_a, r.layer_b, r.layer_c)
            @test layer.outcome == WalkESDTests.LAYER_SKIP
            @test !isempty(layer.reason)
        end
    end

    @test isfile(junit)
    xml = read(junit, String)
    @test occursin("<testsuites>", xml)
    @test occursin("<testsuite name=\"ESD Walker\"", xml)
    @test occursin("tests=\"9\"", xml)
    @test occursin("failures=\"0\"", xml)
    @test occursin("skipped=\"9\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform\"", xml)
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
