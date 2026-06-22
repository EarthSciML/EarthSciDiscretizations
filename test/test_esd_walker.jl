using Test
using TestItems

# Tests for the multi-layer CI walker that validates rule files under
# discretizations/. The walker discovers rules via load_rules and runs five
# layers per rule: (A) canonical-form byte-diff via ESS's `discretize`
# rule engine, (B) MMS convergence via the per-topology canonical-pipeline
# runners (`discretize → ArrayOp → eval`; 1d_cartesian_periodic is
# ArrayOp-native, the other implemented families still scalarize per cell;
# families without a runner SKIP naming the missing piece), (B')
# TVD/monotonicity for slope-ratio limiters, (C) integration benchmarks,
# (D) discrete conservation for finite-volume rules. Layer A skips until
# canonical-form fixtures are authored; Layer C is gated on
# ESD_RUN_INTEGRATION=1 and skipped by default; Layer D skips for rules
# without a conservation/ fixture.

@testitem "walker: discovers seeded rules; per-layer outcomes match the catalog ledger" begin
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
            "centered_2nd_deriv_uniform",
            "laplacian_2nd_uniform_cartesian",
            "nonlinear_laplacian_uniform",
            "mixed_deriv_2nd_uniform",
            "dirichlet_bc",
            "neumann_bc",
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

    # Layer A: passes for the rules that ship a canonical/ fixture
    # (centered_2nd_uniform, centered_2nd_uniform_vertical,
    # centered_2nd_uniform_latlon, nonlinear_laplacian_uniform,
    # mixed_deriv_2nd_uniform, and — since the dsc-kswm follow-ups
    # authored byte contracts through the live pipeline — upwind_1st and
    # centered_2nd_deriv_uniform (scheme + use: form),
    # laplacian_2nd_uniform_cartesian (2D arrayop lift + periodic-fold
    # contract), and weno5_advection (promoted from its applicable:false
    # stub)); skips for every other rule with reason "no canonical or
    # rewrite fixtures" (dsc-aez introduced the rewrite/ variant; no rule
    # has one committed yet — see the synthetic-rule unit tests below).
    # Layer C always skips unless ESD_RUN_INTEGRATION=1.
    # Layer B: per-topology canonical-pipeline runners (dsc-kswm) drive
    # `discretize → ArrayOp → eval` for the implemented families
    # (1d_cartesian_periodic ArrayOp-native; vertical/latlon/arakawa still
    # per-cell-scalar pending ESS scheme dispatch beyond the cartesian
    # foundation). Rules in families without a runner SKIP with the
    # unified `_LAYER_B_PIPELINE_PENDING` reason; rules whose convergence
    # fixture declares `applicable:false` keep their existing
    # `fixture-declared not applicable` SKIP path.
    # The pending-canonical set is EMPTY as of the dsc-kswm follow-ups:
    # every rule with an applicable:true convergence fixture now has a
    # Layer-B runner. History (rules that once sat here): the esd-bbp
    # vertical/latlon/arakawa activations, upwind_1st (esd-0ip),
    # laplacian_2nd_uniform_cartesian (dsc-vst2), and
    # centered_2nd_deriv_uniform (whose second-derivative MMS landed with
    # the pending-set retirement). A rule re-enters by shipping an
    # applicable:true convergence fixture for a topology family without a
    # runner.
    pending_canonical_layer_b = Set(Tuple{String, String}[])
    # vertical_remap (dsc-otd) is structurally a phase-hook operation (Lagrangian
    # → Eulerian re-gridding between timesteps), not a §7 stencil rule. The
    # JSON file is retained as a reference artifact documenting the Lin (2004)
    # PPM-remap math and AST shape for the eventual phase-hook contract; there
    # is no imperative implementation in any binding. The convergence fixture
    # ships applicable:false with a phase-hook deferral reason so the walker
    # SKIPs Layer-B with "fixture-declared not applicable".
    # flux_1d_ppm canonical fixture is activated (applicable:true); Layer-A passes
    # via the inline pattern+replacement ESM document and is handled by the
    # explicit elseif below (Layer-A PASS, Layer-B SKIP).
    not_applicable_layer_b = Set(
        [
            # centered_2nd_nonuniform_vertical: the 1d_vertical_column runner
            # loads the replacement from a canonical/ fixture (absent here) and
            # binds uniform h only — per-cell dz[k] bindings pending dsc-yz0m.
            ("finite_difference", "centered_2nd_nonuniform_vertical"),
            ("finite_difference", "periodic_bc"),
            ("finite_volume", "flux_limiter_minmod"),
            ("finite_volume", "flux_limiter_superbee"),
            ("finite_volume", "lax_friedrichs_flux"),
            ("finite_volume", "vertical_remap"),
        ]
    )
    # Layer B' (limiter): TVD slope-ratio limiters ship a `monotonicity/`
    # fixture kind under discretizations/<rule>/fixtures/ and the walker
    # exercises Sweby-region + 1D advection TVD checks against the rule's
    # AST (dsc-8vu). All other rules SKIP this layer.
    pass_layer_limiter = Set(
        [
            ("finite_volume", "flux_limiter_minmod"),
            ("finite_volume", "flux_limiter_superbee"),
        ]
    )
    # Layer D (conservation): finite-volume rules that ship a
    # `conservation/` fixture (dsc-559). Currently divergence_arakawa_c
    # exercises the 2D periodic divergence telescoping check, and
    # flux_limiter_minmod exercises the 1D MUSCL telescoping check. All
    # other rules SKIP because they have no conservation/ fixture.
    pass_layer_d = Set(
        [
            ("finite_volume", "divergence_arakawa_c"),
            ("finite_volume", "flux_limiter_minmod"),
        ]
    )
    # Canonical-skip-only fixtures declaring `applicable:false` that also carry
    # no convergence fixture, so Layer-A SKIPs via `_fixture_applicable_skip`
    # and Layer-B SKIPs with "no convergence fixtures at...".
    canonical_skip_only_fv3 = Set(
        [
            # dirichlet_bc (esd-0eu): rewrite fixture ships `applicable:false`
            # pending the kind/side discrimination flip (layer-A SKIP) and carries
            # no convergence fixture (layer-B SKIP "no convergence fixtures").
            ("finite_difference", "dirichlet_bc"),
            # neumann_bc (esd-6g4.9 / G10, esd-6k1): PROMOTED out of this set. The
            # ESS fn/dim BC-kind matcher (ess-bps/ess-tox, G8) landed, so the
            # rewrite fixture is applicable:true — bc(neumann, xmin, [u, g]) fires
            # the single side-generic neumann_bc rule (bind_side_spacing binds
            # $h = 1/N). Layer-A PASSes; see the explicit elseif below. The numeric
            # INTEGRATION path now solves too (esd-6k1) via the ess-hjg makearray
            # lowering — see test/test_bc_ic_goldens.jl.
            # robin_bc (esd-m9v, esd-6g4.8/G9, esd-6k1): PROMOTED out of this set.
            # The fn/dim+args firing encoding (the three Robin coefficients ride as
            # trailing args [$u, $a, $b, $g]) + the ess-hjg makearray lowering make
            # the rewrite fixture applicable:true and the INTEGRATION path solve.
            # Layer-A PASSes; see the explicit elseif below. Supersedes the
            # esd-6g4.8 ROBIN_BC_INFEASIBILITY verdict.
            # staggered_1st_uniform (esd-6g4.13): promoted out of this set. The
            # arrayop replacement (EINSUM-4) now drives end-to-end through ESS
            # `discretize` once the rule's `applies_to` carries `dim: $x` (the
            # missing field — not "stencil-schema dispatch" — was what blocked the
            # earlier applicable:false stub). Layer-A PASSes via the canonical byte
            # contract; see the explicit elseif below.
        ]
    )
    # The hot-path FV rules (dsc-ntxo, audit dsc-ztvz / F6) ship
    # canonical-skip-only fixtures declaring `applicable:false`. ESS now
    # dispatches single-output cartesian stencils via `use:`-schemes
    # (esm-j1u), but these reconstruction/flux rules need what ESS still
    # lacks: multi-output emission (stencil objects keyed by output edge)
    # and face-staggered bindings — see each fixture's `skip_reason` for
    # the per-rule contract. Layer-A SKIPs via `_fixture_applicable_skip`.
    # Layer-B keeps its existing `_LAYER_B_PIPELINE_PENDING` SKIP (these
    # rules also ship convergence fixtures with applicable:true). Once ESS
    # gains multi-output emission, these stubs flip to PASS the same way
    # as the FV3 ports.
    canonical_skip_only_hot_path = Set(
        [
            # divergence_arakawa_c removed (esd-bbp): now PASSes Layer-B via
            # canonical pipeline with mms_kind="vec_sincos_2d_periodic" (O(h²)).
            # ppm_reconstruction removed (dsc-a7b2): now PASSes Layer-B via the
            # ArrayOp-native fv_cell_average_1d runner (per-output schemes,
            # mms_kind="sin_2pi_x_cell_average", O(h⁴) edge interpolation);
            # its Layer-A canonical stub stays applicable:false pending
            # document-level multi-output emission.
            # weno5_advection removed: now PASSes Layer-B via the
            # ArrayOp-native 1d_cartesian_periodic runner with the auxiliary
            # velocity field bound from _LAYER_B_MMS_AUX
            # (mms_kind="sin_2pi_x_unit_advection", O(h⁵) Jiang-Shu FD-WENO);
            # its Layer-A canonical stub stays applicable:false.
            ("finite_volume", "weno5_advection_2d"),
        ]
    )
    # esd-6g4.13 (G14): high-order grad (4/6/8) + upwind ported to the vertical
    # and latlon grid families. Each ships a canonical byte contract (Layer-A
    # PASS via `discretize`) and no convergence/ fixture (Layer-B SKIP "no
    # convergence fixtures"), mirroring the staggered_1st_uniform /
    # upwind_1st_nonuniform numeric-UNIT pattern. The vertical rules are the
    # uniform-spacing cartesian central stencils with dz=h (coefficients
    # verified by the identical cartesian centered_4th/6th/8th_uniform Layer-B);
    # the latlon rules carry the spherical metric (1/(R·cos_lat·dlon) for the
    # lon term, 1/(R·dlat) for the lat term) on standard central coefficients.
    pass_layer_a_canonical_only_g14 = Set(
        [
            ("finite_difference", "centered_4th_uniform_vertical"),
            ("finite_difference", "centered_6th_uniform_vertical"),
            ("finite_difference", "centered_8th_uniform_vertical"),
            ("finite_difference", "upwind_1st_vertical"),
            ("finite_difference", "centered_4th_uniform_latlon"),
            ("finite_difference", "centered_6th_uniform_latlon"),
            ("finite_difference", "centered_8th_uniform_latlon"),
            ("finite_difference", "upwind_1st_latlon"),
        ]
    )
    for r in results
        @test r.layer_c.outcome == WalkESDTests.LAYER_SKIP
        @test !isempty(r.layer_c.reason)
        key = (String(r.family), r.name)
        if r.family === :finite_difference && r.name == "centered_2nd_uniform_vertical"
            # centered_2nd_uniform_vertical (vertical) ships a canonical/
            # fixture, so Layer A passes via the ESS rule engine (dsc-cjh).
            # Layer B now passes via the canonical pipeline (esd-bbp):
            # replacement AST added to rule JSON, fixture declares
            # mms_kind="sin_2pi_x_periodic", topology key resolves to
            # 1d_vertical_column, runner measures O(h²).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_2nd_uniform"
            # centered_2nd_uniform (cartesian) ships a canonical/ fixture
            # (dsc-3sg) so Layer-A passes via the ESS rule engine. Layer-B
            # passes via the ArrayOp-native canonical pipeline (dsc-kswm):
            # the fixture declares mms_kind="sin_2pi_x_periodic", topology
            # key resolves to 1d_cartesian_periodic, and the runner rides
            # the rule's replacement AST as a document rule through
            # `discretize(lift_1d_arrayop=true) → ArrayOp → build_evaluator`,
            # measuring O(h²).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_2nd_deriv_uniform"
            # centered_2nd_deriv_uniform: Layer-A passes via its canonical
            # byte contract in scheme + use: form. Layer-B passes via the
            # ArrayOp-native 1d_cartesian_periodic runner on the scheme +
            # use: path, mms_kind="sin_2pi_x_second_derivative" (O(h²)
            # 3-point second derivative).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "weno5_advection"
            # weno5_advection: Layer-A passes via its canonical byte
            # contract (promoted from the applicable:false stub once the
            # replacement-form document path landed; the 39KB expected doc
            # is the regression net for the largest AST in the catalog).
            # Layer-B passes via the ArrayOp-native
            # 1d_cartesian_periodic runner: the rule's replacement AST
            # (pattern div($U * $q, dim=$x)) lowers to canonical components,
            # the auxiliary velocity $U binds to a frozen unit field per
            # _LAYER_B_MMS_AUX, and the Jiang-Shu FD-WENO divergence
            # measures O(h⁵) against d(U·u)/dx — clearing the fixture's
            # expected_min_order of 4.7.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "ppm_reconstruction"
            # ppm_reconstruction (dsc-a7b2): Layer-A passes via its canonical
            # byte contract (promoted from applicable:false once the ESS
            # multi_output_stencil document path landed; 3341-byte canonical
            # doc is the regression net for the multi-output AST). Layer-B
            # esd-agh: lower_stencil_to_scheme retired; fv_cell_average_1d runner
            # removed from _LAYER_B_SUPPORTED_TOPOLOGIES. ppm_reconstruction still
            # carries stencil form (EINSUM-8 territory) so Layer-B SKIPs with the
            # canonical-pipeline-pending reason until ppm_reconstruction is migrated
            # to authored replacement form.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("Layer-B awaits canonical-pipeline replacement", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "laplacian_2nd_uniform_cartesian"
            # laplacian_2nd_uniform_cartesian: Layer-A passes via its
            # canonical byte contract — the first to pin the 2D arrayop lift
            # and periodic ifelse folding. Layer-B passes via the
            # ArrayOp-native 2d_cartesian_periodic runner (dsc-vst2): the
            # stencil lowers to canonical i/j component form (RFC §7.1),
            # the fixture declares mms_kind="sin2pix_sin2piy_periodic", the
            # classifier routes on the stencil's two distinct selector axes,
            # and the runner measures O(h²).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "upwind_1st"
            # upwind_1st: Layer-A passes via its canonical byte contract in
            # scheme + use: form (the §7.2.1 production path). Layer-B
            # passes via the ArrayOp-native canonical pipeline (dsc-kswm): the
            # stencil form lowers to an ESS scheme + use: rule via
            # `lower_stencil_to_scheme`, fixture declares
            # mms_kind="sin_2pi_x_periodic", topology key resolves to
            # 1d_cartesian_periodic, runner measures O(h).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_4th_uniform"
            # centered_4th_uniform: Layer-A passes via its canonical byte contract.
            # Layer-B passes via the 1d_cartesian_periodic runner: fixture declares
            # mms_kind="sin_2pi_x_periodic" on N=16,32,64,128; the 4th-order FD
            # gradient measures O(h⁴).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_6th_uniform"
            # centered_6th_uniform: Layer-A passes via its canonical byte contract.
            # Layer-B passes via the 1d_cartesian_periodic runner: fixture declares
            # mms_kind="sin_2pi_x_periodic" on N=32,64,128,256; the 6th-order FD
            # gradient measures O(h⁶).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_8th_uniform"
            # centered_8th_uniform: Layer-A passes via its canonical byte contract.
            # Layer-B passes via the 1d_cartesian_periodic runner: fixture declares
            # mms_kind="sin_2pi_x_periodic" on N=16,32,64,128; the 8th-order FD
            # gradient measures O(h⁸). N=256 excluded (stencil coefficients amplify
            # roundoff past the truncation crossover near N~160).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_4th_deriv_uniform"
            # centered_4th_deriv_uniform: Layer-A passes via its canonical byte
            # contract. Layer-B passes via the 1d_cartesian_periodic runner: fixture
            # declares mms_kind="sin_2pi_x_second_derivative" on N=16,32,64,128; the
            # 4th-order FD second derivative measures O(h⁴).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_6th_deriv_uniform"
            # centered_6th_deriv_uniform: Layer-A passes via its canonical byte
            # contract. Layer-B passes via the 1d_cartesian_periodic runner: fixture
            # declares mms_kind="sin_2pi_x_second_derivative" on N=16,32,64,128; the
            # 6th-order FD second derivative measures O(h⁶). N=256 excluded (stencil
            # coefficients amplify roundoff past the truncation crossover near N~320,
            # observed as order 3.5 on the N=128→256 step at N=32,64,128,256).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_8th_deriv_uniform"
            # centered_8th_deriv_uniform: Layer-A passes via its canonical byte
            # contract. Layer-B passes via the 1d_cartesian_periodic runner: fixture
            # declares mms_kind="sin_2pi_x_second_derivative" on N=16,32,64 (2
            # convergence steps); the 8th-order FD second derivative measures O(h⁸).
            # N>=128 excluded: large stencil coefficients (max 14350) amplify roundoff
            # as eps/(h²*5040), dominating truncation past the crossover near N~76.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "nonlinear_laplacian_uniform"
            # nonlinear_laplacian_uniform (esd-1p7) ships a canonical/ fixture
            # so Layer-A passes via the ESS rule engine. Layer-B passes via the
            # ArrayOp-native 1d_cartesian_periodic runner: the flux-form
            # pattern grad($f * grad($u, $x), $x) binds the coefficient $f to
            # a frozen space-varying field (f = 2 + sin(2πx)) per
            # _LAYER_B_MMS_AUX, mms_kind="sin_2pi_x_nonlinear_diffusion",
            # and the conservative discretization measures O(h²).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "mixed_deriv_2nd_uniform"
            # mixed_deriv_2nd_uniform (esd-wdv) ships a canonical/ fixture so
            # Layer-A passes via the ESS rule engine. Layer-B passes via the
            # ArrayOp-native 2d_cartesian_periodic runner: the nested pattern
            # grad(grad($u, dim=$y), dim=$x) lowers with both dim pattern
            # variables mapped to canonical components ($x → i, $y → j),
            # mms_kind="sin2pix_sin2piy_mixed_deriv", and the 4-corner cross
            # stencil measures O(h²).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_2nd_uniform_latlon"
            # centered_2nd_uniform_latlon (latlon): the canonical/ fixture was
            # activated in esd-9qs (applicable:false dropped once ess-sra landed
            # the per-cell metric-binding evaluator), so Layer-A passes via the
            # ESS rule engine. Layer-B passes via the canonical pipeline
            # (esd-bbp): the fixture declares mms_kind="Y_2_0_unit_sphere",
            # topology key resolves to 2d_latlon_sphere, and the runner measures
            # O(h²) on the lat axis (Y_{2,0} is lon-independent so only lat
            # signal matters).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "divergence_arakawa_c"
            # divergence_arakawa_c (arakawa) ships a canonical/ fixture with
            # applicable:false (Layer-A skips). Layer-B now passes via the canonical
            # pipeline (esd-bbp): the fixture declares mms_kind="vec_sincos_2d_periodic",
            # topology key resolves to 2d_arakawa_periodic, and the runner measures O(h²)
            # for the divergence of (sin(2πx)cos(2πy), cos(2πx)sin(2πy)).
            # Layer-D (conservation) continues to pass as before.
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "divergence_mpas"
            # divergence_mpas (esd-6g4.1): the flagship MPAS flux-form divergence
            # and the unstructured-operator TEMPLATE for the G-series. Layer-A
            # passes via its canonical byte contract — an inline
            # pattern+replacement rule (div($F, dim=cell) -> reduction arrayop
            # weighted by edge_sign_on_cell * dv_edge / area_cell) on an MPAS grid
            # with an edge-located F state; ESS discretize lowers div(F) and the
            # walker byte-compares to the pinned expected.esm. Layer-B passes via
            # the new unstructured_divergence runner
            # (mms_kind="div_const_field_sphere"): a constant cartesian field's
            # edge-normal flux F_e=V·n̂_e is injected as the edge state, the
            # builtin Voronoi dual ladder (level 3->4->5) is swept, and the
            # cell-centered divergence is measured against -2(V·r̂)/R, clearing
            # the fixture's expected_min_order of 0.9 (O(h) L∞ — the MPAS C-grid
            # divergence is first-order in L∞ on the irregular icosahedral-Voronoi
            # mesh, exactly like nn_diffusion_mpas).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "gradient_mpas"
            # gradient_mpas (esd-6g4.2): MPAS edge-normal gradient (grad φ)_e =
            # (φ[c2]−φ[c1])/dc_edge[e], the TRiSK cell-scalar→edge-normal gradient
            # (Ringler et al. 2010, Eq. 22) — the edge-OUTPUT counterpart of the
            # divergence_mpas template (host primitive cells_on_edge added in
            # src/ode_problem.jl). Layer-A passes via its canonical byte contract:
            # an inline pattern+replacement rule (grad($phi, dim=edge) → an
            # edge-output arrayop whose free output index ranges over EDGES, no
            # reduction) on an MPAS grid with an edge-located g state; ESS
            # discretize lowers grad(phi) and the walker byte-compares to the
            # pinned expected.esm. Layer-B passes via the new unstructured_gradient
            # runner (mms_kind="grad_linear_field_sphere"): the linear cartesian
            # field φ=V·r is injected as the cell state, the builtin Voronoi dual
            # ladder (level 3→4→5) is swept, and the edge-located gradient output
            # is measured against the analytic edge-normal gradient V·n̂_e, clearing
            # the fixture's expected_min_order of 1.7 (O(h²) L∞ — unlike the
            # first-order flux-form divergence/diffusion, the edge-normal gradient
            # is genuinely second order on the irregular dual).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "divergence_duo"
            # divergence_duo (esd-6g4.3): the DUO icosahedral-triangular flux-form
            # divergence, the triangular-primal sibling of divergence_mpas. Layer-A
            # passes via its canonical byte contract — an inline pattern+replacement
            # rule (div($F, dim=cell) -> reduction arrayop weighted by
            # edge_sign_on_face * dc_edge / tri_area) on a DUO grid with an
            # edge-located F state; ESS discretize lowers div(F) and the walker
            # byte-compares to the pinned expected.esm. Layer-B passes via the
            # unstructured_divergence runner's DUO branch
            # (mms_kind="div_const_field_sphere"): a constant cartesian field's
            # edge-normal flux F_e=V·n̂_e (circumcenter-chord normal) is injected as
            # the edge state, the builtin icosahedral ladder (level 2->3->4) is
            # swept, and the cell-centered divergence is measured at each cell's
            # circumcenter against -2(V·r̂)/R, clearing the expected_min_order of
            # 0.9 (O(h) L∞ — the C-grid divergence is first-order in L∞ on the
            # irregular icosahedral-triangular mesh, like nn_diffusion_duo).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "advection_mpas"
            # advection_mpas (esd-6g4.2): flux-form MPAS scalar advection
            # div(u*q) — the divergence_mpas reduction with the edge flux
            # u_e·(q[c1]+q[c2])/2 (q centered-interpolated to the edge). Matches
            # the weno5_advection pattern div($u*$q). Layer-A passes via its
            # canonical byte contract: an inline location-GUARDED pattern (the
            # var_location_is guards pin q=cell_center, u=edge_normal so the
            # commutative-product operand roles survive ESS canonicalization).
            # Layer-B SKIPs by fixture decree (applicable:false): the unlimited
            # centered C-grid advection is L∞-low-order on the non-centroidal
            # icosahedral-Voronoi dual — the edge-interpolation defect plateaus
            # L∞ at the most distorted cells (verified MMS-independent: both a
            # constant-cartesian-field and a divergence-free rotation MMS plateau),
            # so a clean L∞-convergent sweep is deferred. The rule's correctness
            # rests on the Layer-A byte contract plus the shared divergence_mpas
            # convergence proof. See MPAS_FLUX_VERDICT.md. The operand-order mirror
            # advection_mpas_velocity_first (a sibling key in advection_mpas.json,
            # not separately walked) makes the operator robust to either variable
            # naming when both variants are referenced.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "gradient_duo"
            # gradient_duo (esd-6g4.3): the DUO C-grid normal gradient, the
            # discrete adjoint of divergence_duo (edge OUTPUT from a cell-centered
            # input). Layer-A passes via its canonical byte contract — an inline
            # pattern+replacement rule (grad($q, dim=edge) -> a reduction arrayop
            # over the two incident faces s∈{1,2} with the centered-difference sign
            # (2s-3), the free index running over edges) on a DUO grid; ESS
            # discretize lowers grad(q) and the walker byte-compares to the pinned
            # expected.esm. Layer-B passes via the new unstructured_gradient runner
            # (mms_kind="grad_const_field_sphere"): the linear scalar q=a·r̂ is
            # injected at cell circumcenters, the builtin icosahedral ladder (level
            # 2->3->4) is swept, and the edge-normal gradient (q[f_hi]-q[f_lo])/
            # dv_edge is measured against a·n̂_e, clearing the expected_min_order of
            # 0.9 (O(h) L∞ — the two-point edge-normal difference is first-order on
            # the irregular icosahedral mesh). Observed orders [0.99, 0.99].
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "advection_duo"
            # advection_duo (esd-6g4.3): flux-form scalar transport on the DUO
            # mesh — the divergence of the advective flux q·U, composing the
            # divergence_duo stencil with a centered edge reconstruction. Layer-A
            # passes via its canonical byte contract (the div($U*$q, dim=cell)
            # pattern lowers to the cell-output reduction). Layer-B SKIPs: the
            # convergence fixture ships applicable:false — the flux-divergence
            # STRUCTURE is validated (the q≡1 reduction equals divergence_duo at
            # O(h)), but the centered edge reconstruction (q[i]+q[nbr])/2 is
            # sub-first-order in L∞ on the irregular icosahedral mesh, so the
            # spatially-varying-scalar MMS pends a higher-order reconstruction.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "flux_duo"
            # flux_duo (esd-6g4.3): the edge-normal advective flux U·q_e (the
            # un-summed integrand of advection_duo), an edge-OUTPUT reconstruction
            # like gradient_duo. Layer-A passes via its canonical byte contract
            # (the flux($U,$q, dim=edge) pattern lowers to the edge-output
            # reduction over the two incident faces). Layer-B SKIPs: applicable:
            # false — the edge reconstruction shares advection_duo's open item
            # (sub-first-order centered interpolation on the irregular mesh).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        elseif r.family === :finite_volume && r.name == "flux_1d_ppm"
            # flux_1d_ppm: Layer-A passes via its canonical byte contract —
            # the PPM replacement AST (CW84 limiter + Courant-fraction integral,
            # fully expanded) is provided inline in canonical/input.esm as a
            # pattern+replacement rule; ESS discretize applies the substitution
            # and the walker byte-compares to the pinned expected.esm. Layer-B
            # SKIPs because the convergence fixture ships applicable:false (face-
            # staggered Courant binding + ghost_width contract pending in ESS;
            # see flux_1d_ppm.json schema_gaps).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "nn_diffusion_mpas"
            # nn_diffusion_mpas (esd-cal): unstructured_ode runner on MPAS Voronoi
            # mesh. No canonical/ fixture (Layer-A skips). Layer-B passes via the
            # unstructured_ode runner: mms_kind="sin_lon_cos_lat", builtin Voronoi
            # dual mesh ladder (x1.642 → x1.2562 → x1.10242), min_order=1.9.
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no canonical or rewrite fixtures", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "nn_diffusion_duo"
            # nn_diffusion_duo (esd-cal): unstructured_ode runner on DUO icosahedral
            # mesh. No canonical/ fixture (Layer-A skips). Layer-B passes via the
            # unstructured_ode runner: mms_kind="sin_lon_cos_lat", builtin icosahedral
            # mesh ladder (level 2 → level 3 → level 4), min_order=1.9.
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no canonical or rewrite fixtures", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_PASS
            @test occursin("min order", r.layer_b.reason)
        elseif key in pending_canonical_layer_b
            # The remaining `applicable:true` convergence fixtures: Layer-A
            # is unconstrained here (centered_2nd_uniform's canonical fixture
            # is asserted in the special-case above; divergence_arakawa_c had
            # pre-existing canonical drift before esm-4t5 retired Layer-B's
            # evaluator; the others ship no canonical/ fixture). Layer-B
            # SKIPs with the unified retirement reason. The n_fail tally
            # below absorbs any unrelated Layer-A drift dynamically.
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("Layer-B awaits canonical-pipeline replacement", r.layer_b.reason)
        elseif key in not_applicable_layer_b
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_b.reason)
        elseif key in canonical_skip_only_fv3
            # Layer-A SKIPs via the `_fixture_applicable_skip` honoring
            # of `applicable:false` in `canonical/input.esm`. Layer-B
            # SKIPs with "no convergence fixtures at..." because these
            # ports carry no convergence fixture.
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif key in canonical_skip_only_hot_path
            # Layer-A SKIPs via `_fixture_applicable_skip`. Layer-B keeps
            # the unified `_LAYER_B_PIPELINE_PENDING` SKIP reason because
            # these rules carry an applicable:true convergence fixture.
            @test r.layer_a.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("fixture-declared not applicable", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("Layer-B awaits canonical-pipeline replacement", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "upwind_1st_nonuniform"
            # upwind_1st_nonuniform (esd-7h2): migrated to arrayop replacement form.
            # Layer-A passes via inline pattern+replacement canonical byte contract.
            # Layer-B SKIPs: convergence fixture has no mms_kind registered.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test !isempty(r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "centered_2nd_nonuniform_cartesian"
            # centered_2nd_nonuniform_cartesian (esd-7h2): migrated to arrayop replacement form.
            # Layer-A passes via inline pattern+replacement canonical byte contract.
            # Layer-B SKIPs: convergence fixture has no mms_kind registered.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test !isempty(r.layer_b.reason)
        elseif r.family === :ic && r.name == "expression_ic"
            # expression_ic (esd-6g4.12 / G13): ships a canonical/ fixture that
            # drives the ESS ic-arrayop materializer on the structured grid
            # families — latlon (8×8) and arakawa (4×4) carried as two grids +
            # two models in one document. Layer-A passes via the whole-document
            # byte golden (each model's IC lowers to an arrayop with grid-correct
            # ranges and coord_<dim> substitution). Layer-B/C/D SKIP: an IC is an
            # initialization transform, not a spatial-operator MMS / conservation
            # target, so the rule ships no convergence/ or conservation/ fixture.
            # Unstructured families (duo/mpas) are DECLARATIVE-INFEASIBLE here — see
            # discretizations/ic/expression_ic/UNSTRUCTURED_IC_INFEASIBILITY.md.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "staggered_1st_uniform"
            # staggered_1st_uniform (esd-6g4.13): Layer-A passes via its canonical
            # byte contract. The cc_to_face arrayop replacement drives the MAC
            # staggered gradient grad(p)|cell_center → dv/dt|face_x through ESS
            # `discretize` (single-equation document; the face_to_cc direction is
            # the mirror stencil documented in the rule JSON). Layer-B SKIPs: the
            # rule ships no convergence/ fixture (a 1D staggered MMS runner is
            # out of scope here — the catalog rule's numeric proof is the Layer-A
            # byte contract).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "zero_gradient_bc"
            # zero_gradient_bc (esd-6g4.7 / G7): the homogeneous-Neumann
            # (du/dn=0) BC kind, now firing through the ESS rule engine via the
            # fn/dim BC-kind matcher landed in ess-bps (G8). The rewrite fixture
            # lifts a synthetic `bc` node with fn="zero_gradient", dim="xmin"
            # and rewrites it to index(u, 1) — the nearest interior cell, giving
            # a zero one-sided difference at the boundary face. Layer-A PASSes
            # via the rewrite canonical-form byte match; Layer-B SKIPs (a
            # ghost-cell rewrite carries no MMS convergence fixture).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "neumann_bc"
            # neumann_bc (esd-6g4.9 / G10, re-based by esd-6k1): the nonzero-flux
            # generalization of zero_gradient_bc, firing through the same ESS
            # fn/dim BC-kind matcher (ess-bps/ess-tox, G8). The rewrite fixture
            # lifts a `bc` node with fn="neumann", dim="xmin" and the flux value as
            # the second arg, then rewrites to index(u, 0) + h*value — the LOCAL
            # 0-based first interior cell plus the grid-spacing-scaled flux, with
            # $h = 1/N bound by `bind_side_spacing` (N=4 -> 0.25*g). One
            # side-generic rule (dim=$side) serves both sides; the makearray
            # lowering (ess-hjg) re-indexes index(u,0) per side (min: u[1];
            # max: u[N]). Layer-A PASSes via the rewrite canonical-form byte match;
            # Layer-B SKIPs (a ghost-cell rewrite carries no MMS convergence
            # fixture). The numeric INTEGRATION path (build_ode_problem) now solves
            # too (esd-6k1) via the ess-hjg makearray lowering — see
            # test/test_bc_ic_goldens.jl (nonzero-Neumann 1D du).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("rewrite canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "robin_bc"
            # robin_bc (esd-m9v, esd-6g4.8/G9, activated by esd-6k1): the Robin BC
            # kind αu + β∂u/∂n = γ, firing through the ESS fn/dim BC-kind matcher
            # (ess-bps/ess-tox, G8) plus the makearray-region lowering (ess-hjg).
            # The rewrite fixture lifts a `bc` node with fn="robin", dim="xmin" and
            # the three coefficients riding as trailing args [u, a, b, g], then
            # rewrites to the ghost (2·h·γ + (2·β − α·h)·index(u,0)) / (α·h + 2·β),
            # with $h = 1/N bound by `bind_side_spacing` (N=4 -> 0.25). The LOCAL
            # 0-based index(u,0) is re-indexed per side by ess-hjg. Layer-A PASSes
            # via the rewrite canonical-form byte match; Layer-B SKIPs (a ghost-cell
            # rewrite carries no MMS convergence fixture). The numeric INTEGRATION
            # path (build_ode_problem) solves too (esd-6k1) — see
            # test/test_bc_ic_goldens.jl (Robin 1D du). Supersedes the esd-6g4.8
            # ROBIN_BC_INFEASIBILITY verdict.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("rewrite canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif key in pass_layer_a_canonical_only_g14
            # esd-6g4.13 (G14): vertical / latlon high-order grad + upwind.
            # Layer-A passes via the canonical byte contract; Layer-B SKIPs (no
            # convergence/ fixture — see the set's definition comment above).
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :integral && r.name == "midpoint_1d"
            # midpoint_1d (esd-6g4.14): full-domain 1-D midpoint quadrature shipped
            # as a `reduce` sum_product arrayop. Layer-A passes via its canonical
            # byte contract — the `integral` op rewrites to `dx · Σ_j u[j]` with the
            # reduction bound `index(size_x, 1)` tracking the grid size as a
            # const_array (the GAP-A bypass; see discretizations/integral/
            # INTEGRAL_FEASIBILITY.md). Layer-B SKIPs: MMS convergence measures a
            # differential-operator's algebraic order, which is N/A for a quadrature
            # (midpoint is O(h²) in *quadrature* error; the rule ships no convergence
            # fixture, so the runner SKIPs with "no convergence fixtures").
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
        elseif r.family === :finite_difference && r.name == "interface_bc"
            # interface_bc (esd-6g4.7 / G7): the value-continuity coupling BC
            # kind (u(x*) = v(x*)). The rewrite fixture lifts a `bc` node with
            # fn="interface", dim="xmin" and the coupled variable as the second
            # arg, then rewrites to index(u, N) where N is bound from the grid
            # via the `bind_side_dim_size` guard — byte-identical to ESS's own
            # conformance golden interface_bc_lowering (index(u, 4) on a size-4
            # dim). Layer-A PASSes via the rewrite canonical-form byte match;
            # Layer-B SKIPs (ghost-cell rewrite, no convergence fixture). The
            # `flux_match` interface modifier (esm-spec §11.5) is a DECLARATIVE
            # GAP — see discretizations/finite_difference/interface_bc/
            # FLUX_MATCH_DECLARATIVE_GAP.md.
            @test r.layer_a.outcome == WalkESDTests.LAYER_PASS
            @test occursin("canonical-form match", r.layer_a.reason)
            @test r.layer_b.outcome == WalkESDTests.LAYER_SKIP
            @test occursin("no convergence fixtures", r.layer_b.reason)
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
    # Layer B: esd-0ip lands the 1d_cartesian_periodic runner; esd-bbp extends
    # it to 1d_vertical_column, 2d_latlon_sphere, and 2d_arakawa_periodic;
    # dsc-vst2 adds the ArrayOp-native 2d_cartesian_periodic runner; weno5_advection
    # and nonlinear_laplacian_uniform ride the 1d runner via auxiliary-field
    # bindings; centered_2nd_deriv_uniform and mixed_deriv_2nd_uniform
    # land with the pending-set retirement; esd-ecq adds the
    # cubed_sphere_cross_metric runner; esd-cal adds the unstructured_ode
    # runner for nn_diffusion_mpas and nn_diffusion_duo; esd-3d7 activates the
    # higher-order cartesian FD family (centered_4th/6th/8th_uniform and
    # centered_4th/6th/8th_deriv_uniform); esd-agh retires lower_stencil_to_scheme
    # and removes the fv_cell_average_1d runner from supported topologies
    # (ppm_reconstruction reverts to SKIP pending EINSUM-8). Nineteen rules now PASS:
    # centered_2nd_uniform (O(h²)), upwind_1st (O(h)),
    # centered_2nd_uniform_vertical (O(h²)), centered_2nd_uniform_latlon (O(h²)
    # on lat axis), divergence_arakawa_c (O(h²) div test),
    # laplacian_2nd_uniform_cartesian (O(h²) 5-point), weno5_advection (O(h⁵)
    # FD-WENO advective divergence), centered_2nd_deriv_uniform (O(h²)),
    # mixed_deriv_2nd_uniform (O(h²) cross stencil),
    # nonlinear_laplacian_uniform (O(h²) flux form),
    # nn_diffusion_mpas (O(h²) MPAS Voronoi sphere),
    # nn_diffusion_duo (O(h²) DUO icosahedral sphere),
    # centered_4th_uniform (O(h⁴)), centered_6th_uniform (O(h⁶)),
    # centered_8th_uniform (O(h⁸)), centered_4th_deriv_uniform (O(h⁴)),
    # centered_6th_deriv_uniform (O(h⁶)), centered_8th_deriv_uniform (O(h⁸)).
    # All remaining rules with applicable:true convergence fixtures continue to
    # SKIP with `_LAYER_B_PIPELINE_PENDING` pending per-topology follow-up beads.
    layer_b_passes = sum(
        1 for r in results
            if r.layer_b.outcome == WalkESDTests.LAYER_PASS;
        init = 0
    )
    # Two layer-B' (limiter) cases pass (minmod, superbee). All other rules
    # SKIP that layer because they have no monotonicity/ fixture directory.
    layer_limiter_passes = sum(
        1 for r in results
            if (String(r.family), r.name) in pass_layer_limiter;
        init = 0
    )
    # Two layer-D (conservation) cases pass (divergence_arakawa_c,
    # flux_limiter_minmod). All other rules SKIP because no conservation/
    # fixture exists.
    layer_d_passes = sum(
        1 for r in results
            if (String(r.family), r.name) in pass_layer_d; init = 0
    )
    # esd-6g4.1 adds divergence_mpas (O(h) MPAS flux-form divergence) via the
    # new unstructured_divergence Layer-B runner, taking the PASS count to 19.
    # esd-6g4.2 adds gradient_mpas (O(h²) MPAS edge-normal gradient) via the new
    # unstructured_gradient runner, taking the PASS count to 20 (advection_mpas
    # ships Layer-A only — its Layer-B convergence fixture is applicable:false).
    # esd-6g4.3 adds divergence_duo (O(h) DUO icosahedral-triangular flux-form
    # divergence) via the same unstructured_divergence runner's DUO branch,
    # taking the PASS count to 21, and gradient_duo (O(h) DUO C-grid normal
    # gradient) via the unstructured_gradient runner's DUO branch, taking it
    # to 22.
    @test layer_b_passes == 22
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
    @test occursin("classname=\"finite_difference.nonlinear_laplacian_uniform\"", xml)
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

@testitem "walker: BC-kind rules (zero_gradient/interface) fire + discriminate via fn/dim" begin
    # esd-6g4.7 / G7: exercises BOTH sides of the real zero_gradient_bc and
    # interface_bc rule files through the ESS engine, and — the whole point of
    # the fn/dim BC-kind matcher (G8, ess-bps) — proves kind DISCRIMINATION:
    # a `bc` node tagged fn="interface" is NOT rewritten by the zero_gradient
    # rules and vice versa. The walker's per-rule rewrite/ fixture only covers
    # one side per rule; this covers the guarded ($N-binding) and guard-free
    # sides of each, plus the negative cases.
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations
    import EarthSciSerialization
    import JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fd = joinpath(repo_root, "discretizations", "finite_difference")

    # Rewrite `expr_dict` (lifted exactly as a rewrite/ fixture would be, via
    # the fn-aware `_expr_from_json`) with the rules in `rule_path` under
    # `ctx_dict`, returning the canonical JSON of the result.
    rw = function (rule_path, expr_dict, ctx_dict)
        rule_doc = JSON.parse(read(rule_path, String))
        rules = EarthSciSerialization.parse_rules(rule_doc["rules"])
        ctx = WalkESDTests._build_rule_context(ctx_dict)
        e = WalkESDTests._expr_from_json(expr_dict)
        out = EarthSciSerialization.rewrite(e, rules, ctx)
        return EarthSciSerialization.canonical_json(out)
    end

    zg = joinpath(fd, "zero_gradient_bc.json")
    itf = joinpath(fd, "interface_bc.json")
    ctx8 = Dict(
        "grids" => Dict(
            "g1" => Dict(
                "spatial_dims" => ["x"],
                "dim_sizes" => Dict("x" => 8)
            )
        ),
        "variables" => Dict("u" => Dict("grid" => "g1"))
    )
    ctx5 = Dict(
        "grids" => Dict(
            "g1" => Dict(
                "spatial_dims" => ["x"],
                "dim_sizes" => Dict("x" => 5)
            )
        ),
        "variables" => Dict(
            "p" => Dict("grid" => "g1"),
            "q" => Dict("grid" => "g1")
        )
    )

    # zero_gradient: ghost = nearest interior cell. xmin -> index(u,1) (guard
    # free); xmax -> index(u,N) with N bound from the grid (here 8).
    @test rw(
        zg, Dict(
            "op" => "bc", "fn" => "zero_gradient", "dim" => "xmin",
            "args" => ["u"]
        ), ctx8
    ) == "{\"args\":[\"u\",1],\"op\":\"index\"}"
    @test rw(
        zg, Dict(
            "op" => "bc", "fn" => "zero_gradient", "dim" => "xmax",
            "args" => ["u"]
        ), ctx8
    ) == "{\"args\":[\"u\",8],\"op\":\"index\"}"

    # interface: ghost reads the coupled variable's far interior cell. xmax ->
    # index(coupled,1) (guard free); xmin -> index(coupled,N) (N=5 here).
    @test rw(
        itf, Dict(
            "op" => "bc", "fn" => "interface", "dim" => "xmax",
            "args" => ["p", "q"]
        ), ctx5
    ) == "{\"args\":[\"q\",1],\"op\":\"index\"}"
    @test rw(
        itf, Dict(
            "op" => "bc", "fn" => "interface", "dim" => "xmin",
            "args" => ["p", "q"]
        ), ctx5
    ) == "{\"args\":[\"q\",5],\"op\":\"index\"}"

    # Discrimination: the zero_gradient rules must NOT touch an interface bc
    # node, and the interface rules must NOT touch a zero_gradient bc node.
    # An unmatched `bc` node passes through unchanged (still op=="bc").
    itf_node = Dict(
        "op" => "bc", "fn" => "interface", "dim" => "xmin",
        "args" => ["p", "q"]
    )
    zg_node = Dict(
        "op" => "bc", "fn" => "zero_gradient", "dim" => "xmin",
        "args" => ["u"]
    )
    @test occursin("\"op\":\"bc\"", rw(zg, itf_node, ctx5))
    @test occursin("\"op\":\"bc\"", rw(itf, zg_node, ctx8))
    # Sanity: the cross-application really did nothing (no index rewrite).
    @test !occursin("\"op\":\"index\"", rw(zg, itf_node, ctx5))
    @test !occursin("\"op\":\"index\"", rw(itf, zg_node, ctx8))
end

@testitem "walker: neumann_bc nonzero-flux rule fires both sides + discriminates via fn/dim" begin
    # esd-6g4.9 / G10, re-based by esd-6k1: exercises BOTH sides of the real
    # neumann_bc rule file through the ESS engine. The walker's rewrite/ fixture
    # only covers xmin; this also pins xmax. The single side-generic rule
    # (dim=$side, 0-based index(u,0)) fires on both sides — the makearray lowering
    # (ess-hjg) re-indexes per side at integration time — and — the point of the
    # fn/dim matcher (G8) — proves the neumann rule does NOT touch a dirichlet or
    # zero_gradient bc node.
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations
    import EarthSciSerialization
    import JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fd = joinpath(repo_root, "discretizations", "finite_difference")

    rw = function (rule_path, expr_dict, ctx_dict)
        rule_doc = JSON.parse(read(rule_path, String))
        rules = EarthSciSerialization.parse_rules(rule_doc["rules"])
        ctx = WalkESDTests._build_rule_context(ctx_dict)
        e = WalkESDTests._expr_from_json(expr_dict)
        out = EarthSciSerialization.rewrite(e, rules, ctx)
        return EarthSciSerialization.canonical_json(out)
    end

    nm = joinpath(fd, "neumann_bc.json")
    # N=4 grid -> h = 1/N = 0.25; flux carried symbolically as `g`.
    ctx4 = Dict(
        "grids" => Dict(
            "g1" => Dict(
                "spatial_dims" => ["x"],
                "dim_sizes" => Dict("x" => 4)
            )
        ),
        "variables" => Dict("u" => Dict("grid" => "g1"))
    )

    # nonzero-Neumann ghost = local 0-based first interior cell + h*flux. One
    # side-generic rule (dim=$side) emits index(u,0) for BOTH sides; the ess-hjg
    # makearray lowering re-indexes it into the grid frame per side at
    # integration time (min -> u[1]; max -> u[N]). Both carry the grid-spacing-
    # scaled flux 0.25*g (h = 1/4). At value=0 each folds to the homogeneous-
    # Neumann mirror.
    @test rw(
        nm, Dict(
            "op" => "bc", "fn" => "neumann", "dim" => "xmin",
            "args" => ["u", "g"]
        ), ctx4
    ) ==
        "{\"args\":[{\"args\":[\"u\",0],\"op\":\"index\"},{\"args\":[0.25,\"g\"],\"op\":\"*\"}],\"op\":\"+\"}"
    @test rw(
        nm, Dict(
            "op" => "bc", "fn" => "neumann", "dim" => "xmax",
            "args" => ["u", "g"]
        ), ctx4
    ) ==
        "{\"args\":[{\"args\":[\"u\",0],\"op\":\"index\"},{\"args\":[0.25,\"g\"],\"op\":\"*\"}],\"op\":\"+\"}"

    # Discrimination: the neumann rules must NOT touch a dirichlet or a
    # zero_gradient bc node (different fn). Each passes through unchanged.
    dir_node = Dict(
        "op" => "bc", "fn" => "dirichlet", "dim" => "xmin",
        "args" => ["u", "g"]
    )
    zg_node = Dict(
        "op" => "bc", "fn" => "zero_gradient", "dim" => "xmin",
        "args" => ["u"]
    )
    @test occursin("\"op\":\"bc\"", rw(nm, dir_node, ctx4))
    @test occursin("\"op\":\"bc\"", rw(nm, zg_node, ctx4))
    @test !occursin("\"op\":\"index\"", rw(nm, dir_node, ctx4))
    @test !occursin("\"op\":\"index\"", rw(nm, zg_node, ctx4))
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
        @test occursin("advection_1d_convergence", result.reason)
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
        # williamson1_cosine_bell retired to stub in esd-vx3 (cubed_sphere_advection runner removed).
        # All 3 cases are stubs → LAYER_SKIP with "3/3 cases skipped".
        @test result.outcome == WalkESDTests.LAYER_SKIP
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

@testitem "walker: layer C runs 1D diffusion analytic (discrete eigenvalue) via field pipeline (ESD_RUN_INTEGRATION=1)" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    rule = first(filter(r -> r.name == "centered_2nd_deriv_uniform", rules))

    prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    try
        ENV["ESD_RUN_INTEGRATION"] = "1"
        result = WalkESDTests.run_layer_c(rule)
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("diffusion_1d_analytic", result.reason)
        @test occursin("L∞", result.reason)
    finally
        if prior === nothing
            delete!(ENV, "ESD_RUN_INTEGRATION")
        else
            ENV["ESD_RUN_INTEGRATION"] = prior
        end
    end
end

@testitem "walker: layer C runs 2D diffusion analytic (discrete eigenvalue) via field pipeline (ESD_RUN_INTEGRATION=1)" begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    rule = first(filter(r -> r.name == "laplacian_2nd_uniform_cartesian", rules))

    prior = get(ENV, "ESD_RUN_INTEGRATION", nothing)
    try
        ENV["ESD_RUN_INTEGRATION"] = "1"
        result = WalkESDTests.run_layer_c(rule)
        @test result.outcome == WalkESDTests.LAYER_PASS
        @test occursin("diffusion_2d_analytic", result.reason)
        @test occursin("L∞", result.reason)
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
