using Test
using TestItems

# esd-t4h: lower_stencil_to_replacement was deleted. All tests below that used synthetic
# stencil rules or idempotence checks are tombstoned here. Rules that carried authored
# replacements (arakawa, staggered) are now tested directly via Layer-A byte contracts.

# esd-t4h: removed "lower_stencil_to_replacement: passes replacement-form rules through"
# lower_stencil_to_replacement was retired in esd-t4h; synthetic idempotence tests are moot.

# esd-7h2+esd-t4h: removed "lower_stencil_to_replacement: ESS parse_expression accepts lowered upwind_1st AST"
# upwind_1st.json carries an authored replacement; stencil-lowering path retired (esd-t4h).

# esd-7h2: removed "lower_stencil_to_replacement: lowers cartesian upwind_1st_nonuniform (esd-02z)"
# upwind_1st_nonuniform.json was migrated to arrayop replacement form — it no longer has a
# stencil array, so the stencil-lowering path is not exercised by this rule.
# The rule is tested via Layer-A (canonical byte contract) and Layer-B (MMS SKIP, no mms_kind).

# esd-7h2: removed "lower_stencil_to_replacement: ESS parse_expression accepts lowered upwind_1st_nonuniform AST (esd-02z)"
# Same migration as above — rule is now in replacement form, stencil path retired.

# esd-t4h: removed "lower_stencil_to_replacement: errors on missing stencil and replacement"
# esd-t4h: removed "lower_stencil_to_replacement: errors on missing applies_to"
# esd-t4h: removed "lower_stencil_to_replacement: errors on unsupported selector kind"
# lower_stencil_to_replacement was retired; synthetic error-path tests are moot.

# esd-7h2: removed "lower_stencil_to_replacement: errors on axis mismatch"
# The axis-mismatch check was inside _lower_cartesian_entry (retired with the
# cartesian branch). Vertical and latlon kinds validate their own axis fields
# independently; no replacement test needed for the retired cartesian path.

# esd-t4h: removed "lower_stencil_to_replacement: single entry produces bare term (vertical)"
# lower_stencil_to_replacement was retired; synthetic vertical stencil tests are moot.

# esd-t4h: removed "lower_stencil_to_replacement: lowers divergence_arakawa_c"
# esd-t4h: removed "lower_stencil_to_replacement: ESS parse_expression accepts lowered divergence_arakawa_c"
# divergence_arakawa_c carries an authored replacement (esd-eg5); lower_stencil_to_replacement
# was idempotent on it. Rule structure is tested by the Layer-A canonical byte contract.

# esd-t4h: removed "lower_stencil_to_replacement: errors on mixed selector kinds"
# lower_stencil_to_replacement was retired; synthetic mixed-kind error tests are moot.

# esd-7h2: removed "lower_stencil_to_replacement: lowers lax_friedrichs_flux (cartesian)"
# lax_friedrichs_flux.json retains its cartesian stencil (EINSUM-8 scope) but
# _lower_cartesian_entry is retired; the rule is FV and does not flow through
# lower_stencil_to_replacement in any active code path. Tested via Layer-A skip
# (applicable:false in convergence fixture) and the hand-pinned flux test in
# test_lax_friedrichs_flux_rule.jl.

# esd-3d7: removed "lower_stencil_to_replacement: lowers centered_2nd_deriv_uniform (cartesian)"
# centered_2nd_deriv_uniform.json was migrated to arrayop replacement form — it no longer has a
# stencil array, so the stencil-lowering path is not exercised by this rule.
# The rule is tested via Layer-A (canonical byte contract) and Layer-B (MMS convergence).

# esd-t4h: removed "lower_stencil_to_replacement: lowers laplacian_2nd_uniform_cartesian (arakawa)"
# esd-t4h: removed "lower_stencil_to_replacement: lowers staggered_1st_uniform_cc_to_face"
# esd-t4h: removed "lower_stencil_to_replacement: lowers staggered_1st_uniform_face_to_cc"
# These rules carry authored replacements (esd-eg5 / esd-770); lower_stencil_to_replacement
# was idempotent on them. Rule structure is tested by the Layer-A canonical byte contracts.

# =============================================================================
# Latlon family tests
# =============================================================================

# esd-t4h: removed "lower_stencil_to_replacement: lowers centered_2nd_uniform_latlon (latlon)".
# centered_2nd_uniform_latlon now carries an authored replacement AST (no stencil field);
# the catalog contract is exercised by the Layer-A byte contract and Layer-B latlon-sphere runner.

# esd-t4h: removed "lower_stencil_to_replacement: latlon errors on \$-prefixed axis"
# lower_stencil_to_replacement was retired; synthetic latlon stencil tests are moot.

# =============================================================================
# Cubed-sphere family tests
# =============================================================================

@testitem "covariant_laplacian_cubed_sphere: cross_metric schema (esd-p81)" begin
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "covariant_laplacian_cubed_sphere.json",
    )
    raw = JSON.parsefile(path)
    disc = raw["discretizations"]

    # Top-level rule is now a cross_metric composite — not a stencil rule.
    rule = Dict{String, Any}(disc["covariant_laplacian_cubed_sphere"])
    @test rule["kind"] == "cross_metric"
    @test rule["grid_family"] == "cubed_sphere"
    @test !haskey(rule, "stencil")
    @test !haskey(rule, "replacement")

    # 8 cross_metric terms: 2 diagonal second-deriv, 2 cross-deriv, 4 metric-gradient corrections.
    terms = rule["terms"]
    @test length(terms) == 8

    # Term metric_component names match the full covariant Laplacian expansion.
    mc = [t["metric_component"] for t in terms]
    @test "ginv_xi_xi"  in mc
    @test "ginv_eta_eta" in mc
    @test "ginv_xi_eta"  in mc   # two entries (ξη and ηξ for symmetric metric)
    @test "dJgxx_dxi"   in mc
    @test "dJgxe_deta"  in mc
    @test "dJgyy_deta"  in mc
    @test "dJgxe_dxi"   in mc

    # All axis_stencil references resolve to entries in the same discretizations block.
    expected_per_axis = Set([
        "d2_dxi2_cubed_sphere",
        "d2_deta2_cubed_sphere",
        "d2_dxieta_cubed_sphere",
        "d1_dxi_over_J_cubed_sphere",
        "d1_deta_over_J_cubed_sphere",
    ])
    referenced = Set(t["axis_stencil"] for t in terms)
    @test referenced ⊆ expected_per_axis
    for name in referenced
        @test haskey(disc, name)
    end

    # Per-axis stencils are cubed_sphere spec-only; lower_stencil_to_replacement retired (esd-t4h).
    for name in expected_per_axis
        pa = Dict{String, Any}(disc[name])
        @test haskey(pa, "stencil")
        @test !haskey(pa, "replacement")
    end
end

# esd-t4h: removed "lower_stencil_to_replacement: lowers transport_2d (cubed_sphere, finite_volume)".
# lower_stencil_to_replacement was retired in esd-t4h; cubed_sphere stencil rules are spec-only
# and not yet materialised (EINSUM-8 territory). The transport_2d stencil schema is verified by
# the cubed_sphere cross_metric runner tests; lowering is not needed until EINSUM-8 authors replacement.

# =============================================================================
# Vertical family tests
# =============================================================================

# esd-t4h: removed "lower_stencil_to_replacement: lowers centered_2nd_uniform_vertical (vertical)".
# centered_2nd_uniform_vertical now carries an authored replacement AST (no stencil field);
# the catalog contract is exercised by the Layer-A byte contract and Layer-B 1d-vertical-column runner.

# esd-t4h: removed "lower_stencil_to_replacement: vertical errors on bad stagger"
# esd-t4h: removed "lower_stencil_to_replacement: vertical errors on axis mismatch"
# lower_stencil_to_replacement was retired; synthetic vertical stencil error tests are moot.

# esd-3d7: removed "lower_stencil_to_scheme: lowers cartesian upwind_1st to ESS scheme parts" and
# "lower_stencil_to_scheme: ESS parse_schemes + parse_rule accept the lowered parts".
# upwind_1st.json was migrated to arrayop replacement form — it no longer has a stencil array, so
# lower_stencil_to_scheme throws ArgumentError on it. The rule is exercised via Layer-A and Layer-B.

# esd-agh: removed "lower_stencil_to_scheme: errors on non-cartesian grid_family and selector kind"
# lower_stencil_to_scheme is retired in esd-agh; synthetic cartesian stencil error-path tests are moot.

# esd-agh: removed "lower_stencil_to_scheme: lowers unstructured nn_diffusion_mpas to ESS scheme parts"
# nn_diffusion_mpas.json was migrated to authored reduce-arrayop replacement form (esd-agh);
# it no longer has a stencil array, so lower_stencil_to_scheme is not exercised.
# The rule is tested numerically via build_ode_problem tests in test_build_ode_problem.jl.

# esd-agh: removed "lower_stencil_to_scheme: unstructured discretize() produces arrayop equations"
# lower_stencil_to_scheme is retired (esd-agh); the rule now uses authored replacement form directly.

# esd-agh: removed "lower_stencil_to_scheme: multi-output stencil object lowers per output"
# lower_stencil_to_scheme is retired in esd-agh; ppm_reconstruction stencil lowering is moot.
# ppm_reconstruction carries its authored replacement form and is tested via Layer-A byte contract.
