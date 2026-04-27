using Test
using TestItems

# Repo-level tests for the discretizations/ rule catalog. These validate
# that the three canonical rule files (centered_2nd_uniform, upwind_1st,
# periodic_bc) are discoverable by load_rules and carry the expected
# schema-level markers (§7 for schemes, §5.2 for rules; ESS discretization
# RFC). The end-to-end rule-engine exercise (.esm → ESS.parse → ESS.rewrite
# → ESS.evaluate / ESS.verify_mms_convergence) lives in test_esd_walker.jl.

@testitem "centered_2nd_uniform scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 1", content)
end

@testitem "upwind_1st scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "upwind_1st", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # 1st-order upwind uses a one-sided stencil: offsets -1 and 0.
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
end

@testitem "periodic_bc rule is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "periodic_bc", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Periodic BC is a rewrite rule (§5.2), not a scheme (§7): it carries
    # `pattern`/`where`/`replacement` rather than `applies_to`/`stencil`.
    @test occursin("\"pattern\"", content)
    @test occursin("\"where\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"dim_is_periodic\"", content)
    @test occursin("\"mod\"", content)
end

@testitem "centered_2nd_uniform_vertical scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_uniform_vertical", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"vertical\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Vertical centered FD uses per-family selector kind and axis k.
    @test occursin("\"kind\": \"vertical\"", content)
    @test occursin("\"axis\": \"\$k\"", content)
    # Face-staggered MMS dispatch (esm-bhv) — selectors carry the per-face
    # stagger plus an integer offset; the two-point centered stencil reads
    # the cell's own bottom and top faces (offset 0).
    @test occursin("\"stagger\": \"face_bottom\"", content)
    @test occursin("\"stagger\": \"face_top\"", content)
    @test occursin("\"offset\": 0", content)
end

@testitem "centered_2nd_uniform_latlon scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_uniform_latlon", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"latlon\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Latlon centered FD uses literal axis values "lon" / "lat" (per
    # SELECTOR_KINDS.md decision #10) so different stencil entries can carry
    # different metrics. Both axes have offsets -1 and 1.
    @test occursin("\"kind\": \"latlon\"", content)
    @test occursin("\"axis\": \"lon\"", content)
    @test occursin("\"axis\": \"lat\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 1", content)
    # Coefficient symbols: angular spacings dlon/dlat, sphere radius R, and
    # the latitude metric cos_lat (lon-axis only) per decisions #11 and #12.
    @test occursin("\"dlon\"", content)
    @test occursin("\"dlat\"", content)
    @test occursin("\"R\"", content)
    @test occursin("\"cos_lat\"", content)
end

@testitem "covariant_laplacian_cubed_sphere scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "covariant_laplacian_cubed_sphere", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cubed_sphere\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"laplacian\"", content)
    # 9-point covariant stencil: offsets in {-1, 0, 1} along both axes.
    @test occursin("\"axis\": \"xi\"", content)
    @test occursin("\"axis\": \"eta\"", content)
    # 2D in-panel offsets are composed via a `selectors` array (one per axis).
    @test occursin("\"selectors\"", content)
    # Cross-panel ghost handling lives in the grid accessor, NOT the selector:
    # the rule must carry no `panel` field.
    @test !occursin("\"panel\"", content)
    # Metric-tensor bindings the cubed_sphere accessor must supply.
    @test occursin("ginv_xi_xi", content)
    @test occursin("ginv_eta_eta", content)
    @test occursin("ginv_xi_eta", content)
    @test occursin("dJgxx_dxi", content)
    @test occursin("dJgyy_deta", content)
    @test occursin("dJgxe_dxi", content)
    @test occursin("dJgxe_deta", content)
end

@testitem "ppm_edge_cubed_sphere scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_edge_cubed_sphere", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cubed_sphere\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"op\": \"reconstruct_panel_edge\"", content)
    # FV3 two-sided edge formula on the gnomonic equidistant cubed sphere
    # collapses to a 4-cell linear stencil at offsets {-1, 0, 1, 2} along the
    # chosen axis ($x ∈ {xi, eta}); 1D-along-axis selectors composed via a
    # `selectors: [...]` array per SELECTOR_KINDS.md decision #13.
    @test occursin("\"selectors\"", content)
    @test occursin("\"axis\": \"\$x\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 2", content)
    # Cross-panel ghost handling lives in the cubed_sphere grid accessor; the
    # rule must NOT carry a `panel` field.
    @test !occursin("\"panel\"", content)
    # FV3 eq. 6.6 monotonicity clamp is documented separately as a non-linear
    # post-processing step; the rule keeps the linear stencil only.
    @test occursin("monotonicity_constraint", content)
end

@testitem "nn_diffusion_mpas scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "nn_diffusion_mpas", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"unstructured\"", content)
    @test occursin("\"stencil\"", content)
    # Operator class: scalar Laplacian acting at cell centers.
    @test occursin("\"op\": \"laplacian\"", content)
    @test occursin("\"emits_location\": \"cell_center\"", content)
    # Two-row formulation closes ∇²u(c) = Σ_k w_k (u(n_k) − u(c)) — Row 1 is a
    # variable-valence neighbor reduction (decision #6/#7 in SELECTOR_KINDS.md);
    # Row 2 is a self-targeting indirect row whose coeff is the diagonal-weight
    # arrayop sum −Σ_k w_k (decision #6 / RFC §7.2 indirect materialization).
    @test occursin("\"kind\": \"reduction\"", content)
    @test occursin("\"kind\": \"indirect\"", content)
    @test occursin("\"table\": \"cells_on_cell\"", content)
    @test occursin("\"k_bound\": \"k\"", content)
    @test occursin("\"index_expr\": \"\$target\"", content)
    @test occursin("\"arrayop\"", content)
    # Coefficient symbols come from the dsc-7j0 MPAS accessor runtime
    # (SELECTOR_KINDS.md decision #9, snake_case): area_cell, dc_edge, dv_edge,
    # edges_on_cell, n_edges_on_cell.
    @test occursin("\"dv_edge\"", content)
    @test occursin("\"dc_edge\"", content)
    @test occursin("\"area_cell\"", content)
    @test occursin("\"edges_on_cell\"", content)
    @test occursin("\"n_edges_on_cell\"", content)
end

@testitem "rule catalog exposes the seeded finite-difference rules" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    names = Set(r.name for r in rules)
    # Superset assertion: the catalog has grown with grid schemas and other
    # families; only require that the canonical FD rules remain discoverable
    # under :finite_difference.
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
    fd_rules = filter(r -> r.family == :finite_difference, rules)
    fd_names = Set(r.name for r in fd_rules)
    @test "centered_2nd_uniform" in fd_names
    @test "centered_2nd_uniform_vertical" in fd_names
    @test "centered_2nd_uniform_latlon" in fd_names
    @test "nn_diffusion_mpas" in fd_names
    @test "periodic_bc" in fd_names
    @test "upwind_1st" in fd_names
    # finite_volume/ppm_reconstruction (CW84 §1) is the first FV rule.
    @test "ppm_reconstruction" in names
    fv_rules = filter(r -> r.family == :finite_volume, rules)
    @test "ppm_reconstruction" in Set(r.name for r in fv_rules)
    # finite_volume/weno5_advection (Jiang-Shu 1996) — 5th-order WENO flux
    # reconstruction, upwind-biased, 1D uniform Cartesian.
    @test "weno5_advection" in names
    @test "weno5_advection" in Set(r.name for r in fv_rules)
    # FV3 D-grid wind-field rules (dsc-247): vorticity (Stokes' theorem),
    # corner interpolation + Coriolis, kinetic energy, D→C interpolation,
    # and metric-corrected sub-grid sin(α) flux on the cubed sphere.
    for fv3_rule in (
        "fv3_vorticity_cellmean",
        "fv3_vorticity_corner",
        "fv3_absolute_vorticity_cellmean",
        "fv3_kinetic_energy_cell",
        "fv3_d_to_c_xi",
        "fv3_d_to_c_eta",
        "fv3_sinsg_flux_xi",
        "fv3_sinsg_flux_eta",
    )
        @test fv3_rule in names
        @test fv3_rule in Set(r.name for r in fv_rules)
    end
end

@testitem "fv3_vorticity_cellmean scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "fv3_vorticity_cellmean", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"curl\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cubed_sphere\"", content)
    @test occursin("\"stagger\": \"D\"", content)
    @test occursin("\"emits_location\": \"cell_center\"", content)
    @test occursin("\"requires_locations\"", content)
    # D-grid stagger selectors at the four oriented edges of the primal cell
    # (decision #17 in SELECTOR_KINDS.md): u_edge / v_edge with axis xi/eta.
    @test occursin("\"stagger\": \"u_edge\"", content)
    @test occursin("\"stagger\": \"v_edge\"", content)
    @test occursin("\"axis\": \"xi\"", content)
    @test occursin("\"axis\": \"eta\"", content)
    @test occursin("\"offset\": 0", content)
    @test occursin("\"offset\": 1", content)
    # Per-cell metric bindings supplied by the cubed_sphere grid accessor:
    # physical edge lengths and primal-cell area.
    @test occursin("\"dx\"", content)
    @test occursin("\"dy\"", content)
    @test occursin("\"area\"", content)
    # Cubed-sphere stagger enum carries the four staggering symbols.
    @test occursin("\"cubed_sphere_stagger\"", content)
end

@testitem "fv3_vorticity_corner scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "fv3_vorticity_corner", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"emits_location\": \"corner\"", content)
    # 4-point average over the four cells meeting at the corner: in-panel
    # offsets in {-1, 0} along both xi and eta axes (cubed_sphere selectors
    # follow decision #13 — composed via `selectors: [...]` array per entry).
    @test occursin("\"selectors\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
    @test !occursin("\"panel\"", content)  # cross-panel ghost lives in accessor
end

@testitem "fv3_absolute_vorticity_cellmean scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "fv3_absolute_vorticity_cellmean", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume

    content = read(rule.path, String)
    @test occursin("\"op\": \"absolute_vorticity\"", content)
    @test occursin("\"emits_location\": \"cell_center\"", content)
    # ω_abs = ω_rel + 2·Omega_rot·sin(lat) — Coriolis bindings present.
    @test occursin("\"Omega_rot\"", content)
    @test occursin("\"lat\"", content)
    @test occursin("\"op\": \"sin\"", content)
end

@testitem "fv3_kinetic_energy_cell scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "fv3_kinetic_energy_cell", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume

    content = read(rule.path, String)
    @test occursin("\"op\": \"kinetic_energy\"", content)
    @test occursin("\"form\": \"closed_form_bilinear\"", content)
    # Bilinear KE: closed-form `expression` over four named D-grid reads
    # (u_w, u_e, v_s, v_n) — not a sum of linear stencil rows.
    @test occursin("\"reads\"", content)
    @test occursin("\"u_w\"", content)
    @test occursin("\"u_e\"", content)
    @test occursin("\"v_s\"", content)
    @test occursin("\"v_n\"", content)
    @test occursin("\"expression\"", content)
    # Sub-grid metric at cell center (super-grid position 5).
    @test occursin("\"sin_a\"", content)
    @test occursin("\"cos_a\"", content)
end

@testitem "fv3_d_to_c_xi / fv3_d_to_c_eta schemes are discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    for (rname, axis, stagger_out) in (
            ("fv3_d_to_c_xi", "xi", "u_edge"),
            ("fv3_d_to_c_eta", "eta", "v_edge"),
        )
        idx = findfirst(r -> r.name == rname, rules)
        @test idx !== nothing
        rule = rules[idx]
        @test rule.family == :finite_volume

        content = read(rule.path, String)
        @test occursin("\"op\": \"d_to_c\"", content)
        @test occursin("\"dim\": \"$axis\"", content)
        @test occursin("\"emits_location\": \"$stagger_out\"", content)
        # Two-point centered average across the staggered interface.
        @test occursin("\"axis\": \"$axis\"", content)
        @test occursin("\"offset\": -1", content)
        @test occursin("\"offset\": 0", content)
    end
end

@testitem "fv3_sinsg_flux_xi / fv3_sinsg_flux_eta schemes are discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    for (rname, axis, stagger, dl) in (
            ("fv3_sinsg_flux_xi", "xi", "u_edge", "dx"),
            ("fv3_sinsg_flux_eta", "eta", "v_edge", "dy"),
        )
        idx = findfirst(r -> r.name == rname, rules)
        @test idx !== nothing
        rule = rules[idx]
        @test rule.family == :finite_volume

        content = read(rule.path, String)
        @test occursin("\"op\": \"metric_flux\"", content)
        @test occursin("\"dim\": \"$axis\"", content)
        @test occursin("\"form\": \"closed_form_upwind_blended\"", content)
        @test occursin("\"emits_location\": \"$stagger\"", content)
        # Upwind-blended closed form: sin_pos / sin_neg precomputed bindings,
        # |v| via the abs op, multiplied by the physical edge length and dt.
        @test occursin("\"sin_pos\"", content)
        @test occursin("\"sin_neg\"", content)
        @test occursin("\"$dl\"", content)
        @test occursin("\"dt\"", content)
        @test occursin("\"op\": \"abs\"", content)
    end
end
