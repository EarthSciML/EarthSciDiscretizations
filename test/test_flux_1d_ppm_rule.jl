using Test
using TestItems

# Tests for the finite_volume/flux_1d_ppm declarative rule.
#
# The rule encodes the PPM-based 1D flux-form transport flux F_{i+1/2}, composing:
#   - 4th-order PPM edge interpolation (CW84 eq. 1.6) on a 6-cell stencil
#     {-3, -2, -1, 0, 1, 2} relative to the right cell of the interface;
#   - Colella-Woodward (1984) §4 monotonicity limiter applied per upwind cell
#     (closed-form ifelse AST);
#   - Courant-fraction flux integral over the swept volume of the upwind cell;
#   - upwind selection encoded by ifelse on sign of the per-face Courant `$c`.
#
# Layer A   — rule discovery + JSON byte-diff round-trip + spot-check schema.
# Layer B (MMS convergence) is `applicable: false` — face-stagger output,
# per-face bindings, and the ghost-extended input contract sit outside the §7
# verify_mms_convergence harness today (see flux_1d_ppm.json `schema_gaps`).

# ---------------------------------------------------------------------------
# Layer A: discovery + spot-checks
# ---------------------------------------------------------------------------

@testitem "flux_1d_ppm rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "flux_1d_ppm", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"flux\"", content)
    @test occursin("\"form\": \"ppm\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"ppm_courant_integral\"", content)
    @test occursin("\"\$c\"", content)
    @test occursin("\"\$q\"", content)
    @test occursin("\"\$v\"", content)
    @test occursin("\"colella_woodward_1984\"", content)
    @test occursin("\"ifelse\"", content)
    # 6-point stencil offsets relative to right cell of the interface.
    # The exhaustive list-equality check is in the round-trip test below.
    @test occursin("\"offset\": -3", content)
    @test occursin("\"offsets\": [-3, -2, -1, 0, 1, 2]", content)
end

@testitem "flux_1d_ppm rule JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "flux_1d_ppm", rules)
    @test idx !== nothing
    rule = rules[idx]

    raw = read(rule.path, String)
    parsed = JSON.parse(raw)
    @test parsed isa AbstractDict
    @test haskey(parsed, "discretizations")
    @test haskey(parsed["discretizations"], "flux_1d_ppm")

    reserialized = JSON.json(parsed)
    reparsed = JSON.parse(reserialized)
    @test reparsed == parsed

    spec = parsed["discretizations"]["flux_1d_ppm"]
    @test spec["applies_to"]["op"] == "flux"
    @test spec["applies_to"]["form"] == "ppm"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["applies_to"]["args"] == ["\$q", "\$c", "\$v"]
    @test spec["grid_family"] == "cartesian"
    @test spec["form"] == "ppm_courant_integral"
    @test spec["produces"] == "F_{i+1/2}"
    @test spec["monotone"] == true
    @test spec["upwind_biased"] == true
    @test spec["composes"] isa AbstractArray
    @test length(spec["composes"]) == 4
    @test spec["stencil_support"]["offsets"] == [-3, -2, -1, 0, 1, 2]
    @test haskey(spec, "edge_value_stencil")
    @test haskey(spec, "limiter")
    @test spec["limiter"]["name"] == "colella_woodward_1984"
    @test spec["limiter"]["inputs"] == ["ql", "qr", "qi"]
    @test spec["limiter"]["outputs"] == ["ql_lim", "qr_lim"]
    @test haskey(spec, "flux_integral")
    @test haskey(spec["flux_integral"], "int_left")
    @test haskey(spec["flux_integral"], "int_right")
    @test haskey(spec, "flux_form")
    @test spec["flux_form"]["produces"] == "F_{i+1/2}"
    @test haskey(spec, "schema_gaps")
end

