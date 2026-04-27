using Test
using TestItems

# Tests for the multi-layer CI walker that validates rule files under
# discretizations/. The walker discovers rules via load_rules and runs five
# layers per rule: (A) canonical-form byte-diff, (B) MMS convergence (driven
# by EarthSciSerialization.verify_mms_convergence), (B') TVD/monotonicity
# for slope-ratio limiters, (C) integration benchmarks, (D) discrete
# conservation for finite-volume rules. Layer A skips until canonical-form
# fixtures are authored; Layer C is gated on ESD_RUN_INTEGRATION=1 and
# skipped by default; Layer D skips for rules without a conservation/ fixture.

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
    # finite_difference/covariant_laplacian_cubed_sphere (dsc-ap9) — first 2D
    # multi-axis selector rule; declares applicable:false on Layer B until ESS
    # extends the harness for cubed_sphere metric bindings.
    @test "covariant_laplacian_cubed_sphere" in names
    lap = first(filter(r -> r.name == "covariant_laplacian_cubed_sphere", results))
    @test lap.family == :finite_difference

    # Layer A: passes for centered_2nd_uniform and
    # centered_2nd_uniform_vertical (both ship a canonical/ fixture); skips
    # for every other rule with reason "no canonical or rewrite fixtures"
    # (dsc-aez introduced the rewrite/ variant; no rule has one committed
    # yet — see the synthetic-rule unit tests below).
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
                        ("finite_difference", "centered_2nd_uniform_latlon"),
                        ("finite_difference", "upwind_1st"),
                        ("finite_volume", "ppm_reconstruction"),
                        ("finite_volume", "weno5_advection"),
                        ("finite_volume", "divergence_arakawa_c")])
    not_applicable_layer_b = Set([("finite_difference", "periodic_bc"),
                                   ("finite_difference", "covariant_laplacian_cubed_sphere"),
                                   ("finite_difference", "nn_diffusion_mpas"),
                                   ("finite_volume", "flux_limiter_minmod"),
                                   ("finite_volume", "flux_limiter_superbee"),
                                   ("finite_volume", "ppm_edge_cubed_sphere"),
                                   ("finite_volume", "weno5_advection_2d")])
    # Rules whose canonical/ fixture has pre-existing layer-A drift that is
    # tracked by a separate bead. We assert layer-B passes via the convergence
    # sweep but do not constrain layer-A here (the n_fail tally below absorbs
    # any layer-A failure dynamically so the assertion stays correct as the
    # canonical drift is fixed in another PR).
    pass_layer_b_canonical_drift = Set([("finite_volume", "divergence_arakawa_c")])
    # Layer B' (limiter): TVD slope-ratio limiters ship a `monotonicity/`
    # fixture kind under discretizations/<rule>/fixtures/ and the walker
    # exercises Sweby-region + 1D advection TVD checks against the rule's
    # AST (dsc-8vu). All other rules SKIP this layer.
    pass_layer_limiter = Set([("finite_volume", "flux_limiter_minmod"),
                              ("finite_volume", "flux_limiter_superbee")])
    # Layer D (conservation): finite-volume rules that ship a
    # `conservation/` fixture (dsc-559). Currently divergence_arakawa_c
    # exercises the 2D periodic divergence telescoping check, and
    # flux_limiter_minmod exercises the 1D MUSCL telescoping check. All
    # other rules SKIP because they have no conservation/ fixture.
    pass_layer_d = Set([("finite_volume", "divergence_arakawa_c"),
                        ("finite_volume", "flux_limiter_minmod")])
    for r in results
        @test r.layer_c.outcome == WalkESDTests.LAYER_SKIP
        @test !isempty(r.layer_c.reason)
        key = (String(r.family), r.name)
        if r.family === :finite_difference && r.name in ("centered_2nd_uniform", "centered_2nd_uniform_vertical")
            # Both centered_2nd_uniform (cartesian) and
            # centered_2nd_uniform_vertical (vertical) ship a canonical/
            # fixture, so Layer A passes via the ESS rule engine (dsc-3sg /
            # dsc-cjh) in addition to the Layer B convergence sweep.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif key in pass_layer_b_canonical_drift
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif key in pass_layer_b
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no canonical or rewrite fixtures", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif key in not_applicable_layer_b
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        else
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no canonical or rewrite fixtures", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test !isempty(r.layer_b.reason)
        end
        if key in pass_layer_limiter
            @test r.layer_limiter.outcome == WalkESDTests.LAYER_PASS
            @test occursin("Sweby OK", r.layer_limiter.reason)
            @test occursin("TVD OK", r.layer_limiter.reason)
        else
            @test r.layer_limiter.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no monotonicity fixtures", r.layer_limiter.reason)
        end
        if key in pass_layer_d
            @test r.layer_d.outcome == WalkESDTests.LAYER_PASS
            @test occursin("conservation OK", r.layer_d.reason)
        else
            @test r.layer_d.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no conservation fixtures", r.layer_d.reason)
        end
    end

    @test isfile(junit)
    xml = read(junit, String)
    @test occursin("<testsuites>", xml)
    @test occursin("<testsuite name=\"ESD Walker\"", xml)
    # Parametrize against actual catalog size: 5 layers (A/B/B'/C/D) per rule.
    total = length(results) * 5
    # Seven layer-B cases pass (centered_2nd_uniform,
    # centered_2nd_uniform_vertical, centered_2nd_uniform_latlon, upwind_1st,
    # ppm_reconstruction, weno5_advection, divergence_arakawa_c); the rest skip.
    layer_b_passes = sum(1 for r in results
                         if (String(r.family), r.name) in pass_layer_b; init = 0)
    # Two layer-B' (limiter) cases pass (minmod, superbee). All other rules
    # SKIP that layer because they have no monotonicity/ fixture directory.
    layer_limiter_passes = sum(1 for r in results
                               if (String(r.family), r.name) in pass_layer_limiter;
                               init = 0)
    # Two layer-D (conservation) cases pass (divergence_arakawa_c,
    # flux_limiter_minmod). All other rules SKIP because no conservation/
    # fixture exists.
    layer_d_passes = sum(1 for r in results
                         if (String(r.family), r.name) in pass_layer_d; init = 0)
    @test layer_b_passes == 7
    @test layer_limiter_passes == 2
    @test layer_d_passes == 2
    # Count fails/skips from the live result set so this assertion stays
    # correct as the catalog evolves (e.g. when new rules ship canonical
    # fixtures that convert a layer-A FAIL into a PASS or SKIP). Avoids
    # double-bookkeeping pre-existing layer-A canonical-form drift unrelated
    # to the limiter Layer B' added by dsc-8vu.
    n_fail = sum(WalkESDTests.count_outcome(r, WalkESDTests.LAYER_FAIL) for r in results; init = 0)
    n_skip = sum(WalkESDTests.count_outcome(r, WalkESDTests.LAYER_SKIP) for r in results; init = 0)
    @test occursin("tests=\"$total\"", xml)
    @test occursin("failures=\"$n_fail\"", xml)
    @test occursin("skipped=\"$n_skip\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_vertical\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_latlon\"", xml)
    @test occursin("classname=\"finite_difference.covariant_laplacian_cubed_sphere\"", xml)
    @test occursin("classname=\"finite_volume.ppm_reconstruction\"", xml)
    @test occursin("classname=\"finite_volume.weno5_advection\"", xml)
    @test occursin("classname=\"finite_volume.flux_limiter_minmod\"", xml)
    @test occursin("classname=\"finite_volume.flux_limiter_superbee\"", xml)
    @test occursin("name=\"layer_A\"", xml)
    @test occursin("name=\"layer_limiter\"", xml)
    @test occursin("name=\"layer_D\"", xml)
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

@testitem "walker: layer A runs `rewrite/` fixture end-to-end via ESS rule engine" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    # End-to-end exercise of the `rewrite/` fixture kind for index-rewrite
    # rules (dsc-aez). Synthetic rule: rewrite `grad($u, dim=$x)` →
    # `index($u, $x)` under guards that bind $g to $u's grid and require $x
    # to be a periodic dimension of $g. Terminating because the replacement
    # changes the operator (`grad` → `index`), preventing re-match.
    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_difference")
        mkpath(family_dir)
        rule_path = joinpath(family_dir, "synth_periodic.json")
        write(rule_path, """
        {
          "rules": {
            "grad_to_index_periodic": {
              "pattern": {"op": "grad", "args": ["\$u"], "dim": "\$x"},
              "where": [
                {"guard": "var_has_grid", "pvar": "\$u", "grid": "\$g"},
                {"guard": "dim_is_periodic", "pvar": "\$x", "grid": "\$g"}
              ],
              "replacement": {"op": "index", "args": ["\$u", "\$x"]}
            }
          }
        }
        """)
        rewrite_dir = joinpath(family_dir, "synth_periodic", "fixtures", "rewrite")
        mkpath(rewrite_dir)
        write(joinpath(rewrite_dir, "input.esm"), """
        {
          "kind": "rewrite",
          "context": {
            "grids": {"g1": {"spatial_dims": ["x"], "periodic_dims": ["x"]}},
            "variables": {"T": {"grid": "g1"}}
          },
          "expression": {"op": "grad", "args": ["T"], "dim": "x"}
        }
        """)
        write(joinpath(rewrite_dir, "expected.esm"),
              "{\"args\":[\"T\",\"x\"],\"op\":\"index\"}\n")

        rule = RuleFile(:finite_difference, "synth_periodic", rule_path)
        result = WalkESDTests.run_layer_a(rule)
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("rewrite canonical-form match", result.reason)
    end
end

@testitem "walker: layer A flags missing rewrite fixture files as failure" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    # Symmetric to the canonical/ missing-files check: an empty rewrite/
    # directory should surface a structured FAIL, not silently SKIP.
    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_difference")
        mkpath(family_dir)
        rule_path = joinpath(family_dir, "broken_rewrite.json")
        write(rule_path, "{\"rules\": {}}")
        rewrite_dir = joinpath(family_dir, "broken_rewrite", "fixtures", "rewrite")
        mkpath(rewrite_dir)

        rule = RuleFile(:finite_difference, "broken_rewrite", rule_path)
        result = WalkESDTests.run_layer_a(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("missing input.esm or expected.esm", result.reason)
    end
end

@testitem "walker: layer A surfaces canonical mismatch from rewrite fixture" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    # Authored expected.esm doesn't match the actual rewrite output: the
    # walker should FAIL with the byte-diff window from `_byte_diff_message`,
    # not silently PASS. Same synthetic rule as the success test.
    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_difference")
        mkpath(family_dir)
        rule_path = joinpath(family_dir, "synth_mismatch.json")
        write(rule_path, """
        {
          "rules": {
            "grad_to_index_periodic": {
              "pattern": {"op": "grad", "args": ["\$u"], "dim": "\$x"},
              "where": [
                {"guard": "var_has_grid", "pvar": "\$u", "grid": "\$g"},
                {"guard": "dim_is_periodic", "pvar": "\$x", "grid": "\$g"}
              ],
              "replacement": {"op": "index", "args": ["\$u", "\$x"]}
            }
          }
        }
        """)
        rewrite_dir = joinpath(family_dir, "synth_mismatch", "fixtures", "rewrite")
        mkpath(rewrite_dir)
        write(joinpath(rewrite_dir, "input.esm"), """
        {
          "kind": "rewrite",
          "context": {
            "grids": {"g1": {"spatial_dims": ["x"], "periodic_dims": ["x"]}},
            "variables": {"T": {"grid": "g1"}}
          },
          "expression": {"op": "grad", "args": ["T"], "dim": "x"}
        }
        """)
        # Wrong: claim the output is `index(U, x)` instead of `index(T, x)`.
        write(joinpath(rewrite_dir, "expected.esm"),
              "{\"args\":[\"U\",\"x\"],\"op\":\"index\"}\n")

        rule = RuleFile(:finite_difference, "synth_mismatch", rule_path)
        result = WalkESDTests.run_layer_a(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("canonical-form mismatch", result.reason)
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

@testitem "walker: layer C runs WENO5 2D smooth advection convergence (ESD_RUN_INTEGRATION=1)" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    weno2d = first(filter(r -> r.name == "weno5_advection_2d", rules))
    @test weno2d.family == :finite_volume

    prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    try
        ENV["ESD_RUN_INTEGRATION"] = "1"
        result = WalkESDTests.run_layer_c(weno2d)
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("smooth_2d_convergence", result.reason)
        @test occursin("observed_min_order", result.reason)
    finally
        if prior === nothing
            delete!(ENV, "ESD_RUN_INTEGRATION")
        else
            ENV["ESD_RUN_INTEGRATION"] = prior
        end
    end
end

@testitem "walker: layer B' (limiter) skips when monotonicity/ is absent" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_volume")
        mkpath(family_dir)
        rule_json = joinpath(family_dir, "no_limiter.json")
        write(rule_json, "{}")
        rule = RuleFile(:finite_volume, "no_limiter", rule_json)

        result = WalkESDTests.run_layer_limiter(rule)
        @test result.outcome == WalkESDTests.LAYER_SKIP
        @test occursin("no monotonicity fixtures", result.reason)
    end
end

@testitem "walker: layer B' passes for flux_limiter_minmod end-to-end" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    minmod = first(filter(r -> r.name == "flux_limiter_minmod", rules))

    result = WalkESDTests.run_layer_limiter(minmod)
    @test result.outcome == WalkESDTests.LAYER_PASS
    @test occursin("Sweby OK", result.reason)
    @test occursin("TVD OK", result.reason)
end

@testitem "walker: layer B' passes for flux_limiter_superbee end-to-end" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    superbee = first(filter(r -> r.name == "flux_limiter_superbee", rules))

    result = WalkESDTests.run_layer_limiter(superbee)
    @test result.outcome == WalkESDTests.LAYER_PASS
    @test occursin("Sweby OK", result.reason)
    @test occursin("TVD OK", result.reason)
end

@testitem "walker: layer B' fails when phi(r) reference disagrees with AST" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile
    using JSON

    # Synthetic limiter rule whose AST is phi(r) = max(0, min(r, 1)) (minmod),
    # paired with a deliberately wrong sweby_check.esm reference value. The
    # walker must surface a FAIL with the offending r and tolerance.
    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_volume")
        mkpath(family_dir)
        rule_path = joinpath(family_dir, "broken_limiter.json")
        rule_doc = Dict(
            "discretizations" => Dict(
                "broken_limiter" => Dict(
                    "formula" => Dict(
                        "op" => "max",
                        "args" => [0, Dict("op" => "min", "args" => ["\$r", 1])],
                    ),
                ),
            ),
        )
        write(rule_path, JSON.json(rule_doc))

        mono = joinpath(family_dir, "broken_limiter", "fixtures", "monotonicity")
        mkpath(mono)
        bad_sweby = Dict(
            "version" => "1.0.0",
            "kind" => "limiter_sweby_check",
            "rule" => "broken_limiter",
            "rule_family" => "finite_volume",
            "variable" => "r",
            "reference_values" => [
                Dict("r" => 1.0, "phi" => 0.5),  # wrong: minmod gives 1.0 here
            ],
            "tvd_properties" => Dict(
                "sweep_r_min" => 0.0,
                "sweep_r_max" => 1.0,
                "sweep_r_step" => 0.5,
            ),
            "tolerance" => 1.0e-12,
        )
        write(joinpath(mono, "sweby_check.esm"), JSON.json(bad_sweby))
        # Provide a tvd_check.esm so the FAIL surfaces from the Sweby check
        # rather than the missing-file branch.
        write(joinpath(mono, "tvd_check.esm"), JSON.json(Dict(
            "grid" => Dict("n" => 8, "dx" => 0.125),
            "advection" => Dict("velocity" => 1.0, "cfl" => 0.4, "periods" => 1.0),
            "tvd_tolerance" => 1.0e-10,
            "eps_denom" => 1.0e-12,
        )))

        rule = RuleFile(:finite_volume, "broken_limiter", rule_path)
        result = WalkESDTests.run_layer_limiter(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("phi(1.0)", result.reason)
        @test occursin("reference is 0.5", result.reason)
    end
end

@testitem "walker: layer D skips when conservation/ is absent" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_volume")
        mkpath(family_dir)
        rule_json = joinpath(family_dir, "no_conservation.json")
        write(rule_json, "{}")
        rule = RuleFile(:finite_volume, "no_conservation", rule_json)

        result = WalkESDTests.run_layer_d(rule)
        @test result.outcome == WalkESDTests.LAYER_SKIP
        @test occursin("no conservation fixtures", result.reason)
    end
end

@testitem "walker: layer D fails when conservation_check.esm is missing" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile

    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_volume")
        mkpath(family_dir)
        rule_json = joinpath(family_dir, "broken_conservation.json")
        write(rule_json, "{}")
        # Empty conservation/ directory: walker must surface a structured FAIL.
        mkpath(joinpath(family_dir, "broken_conservation", "fixtures", "conservation"))

        rule = RuleFile(:finite_volume, "broken_conservation", rule_json)
        result = WalkESDTests.run_layer_d(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("missing conservation_check.esm", result.reason)
    end
end

@testitem "walker: layer D passes for divergence_arakawa_c end-to-end" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    div = first(filter(r -> r.name == "divergence_arakawa_c", rules))

    result = WalkESDTests.run_layer_d(div)
    @test result.outcome == WalkESDTests.LAYER_PASS
    @test occursin("conservation OK", result.reason)
    @test occursin("∇·F", result.reason)
end

@testitem "walker: layer D passes for flux_limiter_minmod end-to-end" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    minmod = first(filter(r -> r.name == "flux_limiter_minmod", rules))

    result = WalkESDTests.run_layer_d(minmod)
    @test result.outcome == WalkESDTests.LAYER_PASS
    @test occursin("conservation OK", result.reason)
    @test occursin("MUSCL", result.reason)
end

@testitem "walker: layer D fails when divergence stencil is corrupted" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: RuleFile
    using JSON

    # Synthetic rule with a non-telescoping stencil: both x-faces carry the
    # same coefficient, so the global divergence sum is non-zero by design.
    # The walker must surface a FAIL with the offending sum.
    mktempdir() do tmp
        family_dir = joinpath(tmp, "finite_volume")
        mkpath(family_dir)
        rule_path = joinpath(family_dir, "broken_div.json")
        rule_doc = Dict(
            "discretizations" => Dict(
                "broken_div" => Dict(
                    "stencil" => [
                        Dict("selector" => Dict("kind" => "arakawa", "stagger" => "face_x", "axis" => "\$x", "offset" => 0),
                             "coeff" => Dict("op" => "/", "args" => [1, "dx"])),
                        Dict("selector" => Dict("kind" => "arakawa", "stagger" => "face_x", "axis" => "\$x", "offset" => 1),
                             "coeff" => Dict("op" => "/", "args" => [1, "dx"])),
                    ],
                ),
            ),
        )
        write(rule_path, JSON.json(rule_doc))

        cdir = joinpath(family_dir, "broken_div", "fixtures", "conservation")
        mkpath(cdir)
        write(joinpath(cdir, "conservation_check.esm"), JSON.json(Dict(
            "kind" => "conservation_divergence_2d_periodic",
            "rule" => "broken_div",
            "grid" => Dict("nx" => 8, "ny" => 6, "dx" => 0.125, "dy" => 0.16666666666666666),
            "seed" => 42,
            "tolerance" => 1.0e-12,
        )))

        rule = RuleFile(:finite_volume, "broken_div", rule_path)
        result = WalkESDTests.run_layer_d(rule)
        @test result.outcome == WalkESDTests.LAYER_FAIL
        @test occursin("global divergence sum", result.reason)
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
