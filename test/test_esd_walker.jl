using Test
using TestItems

# Tests for the multi-layer CI walker that validates rule files under
# discretizations/. The walker discovers rules via load_rules and runs five
# layers per rule: (A) canonical-form byte-diff via ESS's `discretize`
# rule engine, (B) MMS convergence (currently SKIPped pending the
# canonical-pipeline replacement of the retired ESS `verify_mms_convergence`
# — esm-4t5 / 2026-04-29 single-pathway directive), (B') TVD/monotonicity
# for slope-ratio limiters, (C) integration benchmarks, (D) discrete
# conservation for finite-volume rules. Layer A skips until canonical-form
# fixtures are authored; Layer C is gated on ESD_RUN_INTEGRATION=1 and
# skipped by default; Layer D skips for rules without a conservation/ fixture.

@testitem "walker: discovers seeded rules; layer B SKIPs pending canonical-pipeline replacement" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations
    using JSON

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
            "centered_2nd_deriv_uniform",
            "laplacian_2nd_uniform_cartesian",
            "nonlinear_laplacian_uniform",
            "dirichlet_bc",
            "neumann_bc",
            "spherical_laplacian_uniform",
            "robin_bc",
            "centered_2nd_nonuniform_cartesian",
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
    # multi-axis selector rule.
    @test "covariant_laplacian_cubed_sphere" in names
    lap = first(filter(r -> r.name == "covariant_laplacian_cubed_sphere", results))
    @test lap.family == :finite_difference
    # finite_volume/transport_2d (dsc-hdk) — first 2D cubed-sphere FV rule.
    @test "transport_2d" in names
    tr2d = first(filter(r -> r.name == "transport_2d", results))
    @test tr2d.family == :finite_volume

    # Per-rule assertions derived from fixture files — no hardcoded rule-name
    # lists.  Adding a new rule JSON (with or without fixtures) requires ZERO
    # edits to this file: the expected outcome for each layer is inferred from
    # whether the corresponding fixture directory exists and, for layers A/B,
    # whether the fixture declares `applicable: false`.
    for r in results
        # Layer C always skips unless ESD_RUN_INTEGRATION=1.
        @test r.layer_c.outcome == WalkESDTests.LAYER_SKIP
        @test !isempty(r.layer_c.reason)

        fdir = joinpath(dirname(r.path), r.name, "fixtures")
        has_canonical    = isdir(joinpath(fdir, "canonical"))
        has_rewrite      = isdir(joinpath(fdir, "rewrite"))
        has_convergence  = isdir(joinpath(fdir, "convergence"))
        has_monotonicity = isdir(joinpath(fdir, "monotonicity"))
        has_conservation = isdir(joinpath(fdir, "conservation"))

        # ---- Layer A (canonical-form round-trip) ----------------------------
        if has_canonical || has_rewrite
            # Fixture present: layer A must NOT skip with "no canonical or
            # rewrite fixtures" — it either ran (PASS/FAIL) or skipped due
            # to an applicable:false declaration.
            @test !occursin("no canonical or rewrite fixtures", r.layer_a.reason)
        else
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no canonical or rewrite fixtures", r.layer_a.reason)
        end

        # ---- Layer B (MMS convergence) -------------------------------------
        # Layer B wholly SKIPs post-esm-4t5: ESS retired verify_mms_convergence.
        # The replacement (discretize → ArrayOp → official ESS simulation runner)
        # is tracked at ESD/dsc-kswm.
        if has_convergence
            conv_input = joinpath(fdir, "convergence", "input.esm")
            if isfile(conv_input)
                conv_doc = try
                    JSON.parsefile(conv_input)
                catch
                    Dict{String, Any}()
                end
                @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
                if get(conv_doc, "applicable", true) === false
                    @test occursin("fixture-declared not applicable", r.layer_b.reason)
                else
                    @test occursin("Layer-B awaits canonical-pipeline replacement", r.layer_b.reason)
                end
            else
                # convergence/ dir present but missing input.esm → FAIL
                @test r.layer_b.outcome == WalkESDTests.LAYER_FAIL
                @test occursin("missing input.esm", r.layer_b.reason)
            end
        else
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        end

        # ---- Layer B' (limiter monotonicity) --------------------------------
        if has_monotonicity
            @test r.layer_limiter.outcome == WalkESDTests.LAYER_PASS
            @test occursin("Sweby OK", r.layer_limiter.reason)
            @test occursin("TVD OK", r.layer_limiter.reason)
        else
            @test r.layer_limiter.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no monotonicity fixtures", r.layer_limiter.reason)
        end

        # ---- Layer D (discrete conservation) --------------------------------
        if has_conservation
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
    # Layer B is wholly SKIPped post-esm-4t5: ESS retired verify_mms_convergence.
    layer_b_passes = sum(
        1 for r in results
            if r.layer_b.outcome == WalkESDTests.LAYER_PASS;
        init = 0
    )
    n_fail = sum(WalkESDTests.count_outcome(r, WalkESDTests.LAYER_FAIL) for r in results; init = 0)
    n_skip = sum(WalkESDTests.count_outcome(r, WalkESDTests.LAYER_SKIP) for r in results; init = 0)
    @test layer_b_passes == 0
    @test occursin("tests=\"$total\"", xml)
    @test occursin("failures=\"$n_fail\"", xml)
    @test occursin("skipped=\"$n_skip\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_vertical\"", xml)
    @test occursin("classname=\"finite_difference.centered_2nd_uniform_latlon\"", xml)
    @test occursin("classname=\"finite_difference.nonlinear_laplacian_uniform\"", xml)
    @test occursin("classname=\"finite_difference.covariant_laplacian_cubed_sphere\"", xml)
    @test occursin("classname=\"finite_volume.ppm_reconstruction\"", xml)
    @test occursin("classname=\"finite_volume.weno5_advection\"", xml)
    @test occursin("classname=\"finite_volume.flux_limiter_minmod\"", xml)
    @test occursin("classname=\"finite_volume.flux_limiter_superbee\"", xml)
    @test occursin("classname=\"finite_volume.transport_2d\"", xml)
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
        @test occursin("missing input.esm", result.reason)
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
        write(
            rule_path, """
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
            """
        )
        rewrite_dir = joinpath(family_dir, "synth_periodic", "fixtures", "rewrite")
        mkpath(rewrite_dir)
        write(
            joinpath(rewrite_dir, "input.esm"), """
            {
              "kind": "rewrite",
              "context": {
                "grids": {"g1": {"spatial_dims": ["x"], "periodic_dims": ["x"]}},
                "variables": {"T": {"grid": "g1"}}
              },
              "expression": {"op": "grad", "args": ["T"], "dim": "x"}
            }
            """
        )
        write(
            joinpath(rewrite_dir, "expected.esm"),
            "{\"args\":[\"T\",\"x\"],\"op\":\"index\"}\n"
        )

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
        @test occursin("missing input.esm", result.reason)
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
        write(
            rule_path, """
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
            """
        )
        rewrite_dir = joinpath(family_dir, "synth_mismatch", "fixtures", "rewrite")
        mkpath(rewrite_dir)
        write(
            joinpath(rewrite_dir, "input.esm"), """
            {
              "kind": "rewrite",
              "context": {
                "grids": {"g1": {"spatial_dims": ["x"], "periodic_dims": ["x"]}},
                "variables": {"T": {"grid": "g1"}}
              },
              "expression": {"op": "grad", "args": ["T"], "dim": "x"}
            }
            """
        )
        # Wrong: claim the output is `index(U, x)` instead of `index(T, x)`.
        write(
            joinpath(rewrite_dir, "expected.esm"),
            "{\"args\":[\"U\",\"x\"],\"op\":\"index\"}\n"
        )

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
        write(
            joinpath(mono, "tvd_check.esm"), JSON.json(
                Dict(
                    "grid" => Dict("n" => 8, "dx" => 0.125),
                    "advection" => Dict("velocity" => 1.0, "cfl" => 0.4, "periods" => 1.0),
                    "tvd_tolerance" => 1.0e-10,
                    "eps_denom" => 1.0e-12,
                )
            )
        )

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
                        Dict(
                            "selector" => Dict("kind" => "arakawa", "stagger" => "face_x", "axis" => "\$x", "offset" => 0),
                            "coeff" => Dict("op" => "/", "args" => [1, "dx"])
                        ),
                        Dict(
                            "selector" => Dict("kind" => "arakawa", "stagger" => "face_x", "axis" => "\$x", "offset" => 1),
                            "coeff" => Dict("op" => "/", "args" => [1, "dx"])
                        ),
                    ],
                ),
            ),
        )
        write(rule_path, JSON.json(rule_doc))

        cdir = joinpath(family_dir, "broken_div", "fixtures", "conservation")
        mkpath(cdir)
        write(
            joinpath(cdir, "conservation_check.esm"), JSON.json(
                Dict(
                    "kind" => "conservation_divergence_2d_periodic",
                    "rule" => "broken_div",
                    "grid" => Dict("nx" => 8, "ny" => 6, "dx" => 0.125, "dy" => 0.16666666666666666),
                    "seed" => 42,
                    "tolerance" => 1.0e-12,
                )
            )
        )

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
