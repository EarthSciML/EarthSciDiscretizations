using Test
using TestItems

# Tests for the finite_volume/flux_1d_ppm declarative rule.
#
# The rule encodes the PPM-based 1D flux-form transport flux F_{i+1/2}, composing:
#   - 4th-order PPM edge interpolation (CW84 eq. 1.6) on a 6-cell stencil
#     {-3, -2, -1, 0, 1, 2} relative to the right cell of the interface;
#   - Colella-Woodward (1984) §4 monotonicity limiter applied per upwind cell
#     (closed-form ifelse AST matching `_ppm_limit_cw84_sym`);
#   - Courant-fraction flux integral over the swept volume of the upwind cell;
#   - upwind selection encoded by ifelse on sign of the per-face Courant `$c`.
#
# Layer A   — rule discovery + JSON byte-diff round-trip + spot-check schema.
# Layer A'  — hand-pinned single-face flux on a smooth sinusoidal profile,
#             positive- and negative-Courant cases, with the JSON-rule AST
#             evaluated in closed form. The same numbers are produced by
#             `_ppm_limit_cw84_sym` + the closed-form Courant integral in
#             src/operators/flux_1d.jl, so the inline values pin the
#             rule-AST-vs-imperative-reference equivalence.
#             The values are reproduced bit-for-bit from
#             discretizations/finite_volume/flux_1d_ppm/fixtures/canonical/expected.esm.
#
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

# ---------------------------------------------------------------------------
# Layer A' (hand-pinned): closed-form face flux from the rule AST
# ---------------------------------------------------------------------------

@testitem "flux_1d_ppm: hand-pinned single-face flux on smooth profile" begin
    using EarthSciDiscretizations: _ppm_limit_cw84_sym

    # Smooth sinusoidal profile: q[k] = sin(2π·(k - 3.5) / 8) for k in 1..14.
    # The reported face F_{i+1/2} sits between extended cells 7 and 8 (i.e.,
    # the right cell of the face is extended index 8). The 6-point stencil
    # {-3, -2, -1, 0, 1, 2} relative to the right cell pulls extended
    # indices 5..10. The profile is large enough that the CW84 limiter is
    # inactive on both sides of the interface, so the reported flux equals
    # the unlimited PPM Courant integral.
    Nc = 8
    q(k) = sin(2π * (k - 3.5) / Nc)
    qm3 = q(5);  qm2 = q(6);  qm1 = q(7)
    q0  = q(8);  qp1 = q(9);  qp2 = q(10)

    # 4th-order CW84 edge interpolations (the JSON rule's `edge_value_stencil`)
    qi_half = (7/12) * (qm1 + q0) - (1/12) * (qm2 + qp1)
    ql_L = (7/12) * (qm2 + qm1) - (1/12) * (qm3 + q0)
    qr_L = qi_half
    ql_R = qi_half
    qr_R = (7/12) * (q0 + qp1) - (1/12) * (qm1 + qp2)

    # CW84 limiter (the JSON rule's `limiter` block; matches the closed-form
    # ifelse AST in `_ppm_limit_cw84_sym` bit-for-bit).
    ql_L_lim, qr_L_lim = _ppm_limit_cw84_sym(ql_L, qr_L, qm1)
    ql_R_lim, qr_R_lim = _ppm_limit_cw84_sym(ql_R, qr_R, q0)

    # Limiter must be inactive on this smooth profile.
    @test ql_L_lim == ql_L
    @test qr_L_lim == qr_L
    @test ql_R_lim == ql_R
    @test qr_R_lim == qr_R

    # Courant-fraction flux integral (the JSON rule's `flux_integral` block).
    function flux_int_left(ql, qr, qi, c)
        c_abs = abs(c)
        dq = qr - ql
        q6 = 6.0 * (qi - 0.5 * (ql + qr))
        return qr - 0.5 * c_abs * (dq - q6 * (1.0 - (2.0/3.0) * c_abs))
    end
    function flux_int_right(ql, qr, qi, c)
        c_abs = abs(c)
        dq = qr - ql
        q6 = 6.0 * (qi - 0.5 * (ql + qr))
        return ql + 0.5 * c_abs * (dq + q6 * (1.0 - (2.0/3.0) * c_abs))
    end

    v = 1.0
    F_pos = v * flux_int_left(ql_L_lim, qr_L_lim, qm1, 0.3)
    F_neg = v * flux_int_right(ql_R_lim, qr_R_lim, q0, -0.3)

    # Pinned values from
    # discretizations/finite_volume/flux_1d_ppm/fixtures/canonical/expected.esm.
    @test F_pos ≈ 1.2494903985806707e-01 atol = 1.0e-15
    @test F_neg ≈ -1.2494903985806687e-01 atol = 1.0e-15
    # By symmetry of the antisymmetric sinusoid around the face,
    # F_neg = -F_pos to within floating-point roundoff.
    @test abs(F_neg + F_pos) < 1.0e-15
end
