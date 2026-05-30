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

    # dsc-rar: canonical linear-scheme exemplar — the rule lowers grad(u, dim=x)
    # via a closed arrayop replacement in §4.2 ops (no scheme-specific
    # `stencil`/`selector`/`offset` blobs). The replacement is
    #   (u[$x+1] − u[$x−1]) / (2·dx)
    # wrapped in an arrayop over output index $x; BC handling is delegated to
    # downstream BC rules (e.g. periodic_bc) at lowering time.
    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"index\"", content)
    # The closed-form expression is built from §4.2 ops only — no stencil
    # coefficient blobs and no off-spec selector/offset fields.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"selector\"", content)
    @test !occursin("\"offset\"", content)
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

@testitem "upwind_1st_nonuniform scheme is discoverable and well-formed (esd-02z)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "upwind_1st_nonuniform", rules)
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
    # Same one-sided stencil shape as upwind_1st: offsets -1 and 0.
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
    # Non-uniform variant: coefficients index dx per-cell via index("dx", "$x").
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"dx\"", content)
    @test occursin("\"\$x\"", content)
    # No flat scalar "dx" coefficient — must always be indexed.
    @test !occursin("\"args\": [-1, \"dx\"]", content)
    @test !occursin("\"args\": [1, \"dx\"]", content)
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

@testitem "dirichlet_bc rule is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "dirichlet_bc", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Dirichlet BC is a rewrite rule (§5.2 / §9.2): carries `pattern` with
    # `kind:"dirichlet"` and `side:"xmin"` (spec-compliant annotations; ESS
    # kind/side pattern matching is pending), `replacement` giving the ghost-cell
    # value 2*u_bc - u[0], and `produces` declaring the ghost_var per §9.4.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"dirichlet\"", content)
    @test occursin("\"side\": \"xmin\"", content)
    @test occursin("\"op\": \"bc\"", content)
    # Ghost-cell formula: u_ghost = 2*u_bc - u[0].
    @test occursin("\"op\": \"-\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"op\": \"index\"", content)
    # §9.4 ghost_var produces declaration.
    @test occursin("\"produces\"", content)
    @test occursin("\"ghost_var\"", content)
    # Ghost variable named per §9.4 scheme__logical__side convention.
    @test occursin("dirichlet_bc__", content)
    @test occursin("__xmin", content)
end

@testitem "neumann_bc rule is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "neumann_bc", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Neumann BC is a rewrite rule (§5.2 / §9.2): carries `pattern` with
    # `kind:"neumann"` and `side:"xmin"` (spec-compliant annotations; ESS
    # kind/side pattern matching is pending), `replacement` giving the ghost-cell
    # value u[0] + dx*value (where value is du/dn at xmin, the outward normal
    # derivative), and `produces` declaring the ghost_var per §9.4.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"neumann\"", content)
    @test occursin("\"side\": \"xmin\"", content)
    @test occursin("\"op\": \"bc\"", content)
    # Ghost-cell formula: u_ghost = u[0] + dx*value.
    @test occursin("\"op\": \"+\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"dx\"", content)
    # §9.4 ghost_var produces declaration.
    @test occursin("\"produces\"", content)
    @test occursin("\"ghost_var\"", content)
    # Ghost variable named per §9.4 scheme__logical__side convention.
    @test occursin("neumann_bc__", content)
    @test occursin("__xmin", content)
end

@testitem "robin_bc rule is discoverable and well-formed (esd-m9v)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "robin_bc", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Robin BC is a rewrite rule (§5.2 / §9.2): carries `pattern` with
    # `kind:"robin"` and `side:"xmin"`, `robin_alpha/beta/gamma` coefficients
    # for αu + β∂u/∂n = γ, `replacement` giving the ghost-cell value
    # u_ghost = (2·dx·γ + (2·β - α·dx)·u[0]) / (α·dx + 2·β),
    # and `produces` declaring the ghost_var per §9.4.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"robin\"", content)
    @test occursin("\"side\": \"xmin\"", content)
    @test occursin("\"op\": \"bc\"", content)
    # Coefficients: robin_alpha (α), robin_beta (β), robin_gamma (γ).
    @test occursin("\"robin_alpha\"", content)
    @test occursin("\"robin_beta\"", content)
    @test occursin("\"robin_gamma\"", content)
    # Ghost-cell formula uses division, addition, subtraction, multiplication,
    # index into u[0], and the grid spacing dx.
    @test occursin("\"op\": \"/\"", content)
    @test occursin("\"op\": \"+\"", content)
    @test occursin("\"op\": \"-\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"dx\"", content)
    # §9.4 ghost_var produces declaration.
    @test occursin("\"produces\"", content)
    @test occursin("\"ghost_var\"", content)
    # Ghost variable named per §9.4 scheme__logical__side convention.
    @test occursin("robin_bc__", content)
    @test occursin("__xmin", content)
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
            "upwind_1st_nonuniform",
            "dirichlet_bc",
            "neumann_bc",
            "robin_bc",
            "mixed_deriv_2nd_nonuniform",
            "nonlinear_laplacian_nonuniform",
            "spherical_laplacian_nonuniform",
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
    @test "upwind_1st_nonuniform" in fd_names
    @test "dirichlet_bc" in fd_names
    @test "neumann_bc" in fd_names
    @test "robin_bc" in fd_names
    @test "mixed_deriv_2nd_nonuniform" in fd_names
    @test "nonlinear_laplacian_nonuniform" in fd_names
    @test "spherical_laplacian_nonuniform" in fd_names
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
    # finite_volume/fv3_lin_rood_advection (esd-f4i) — Phase 3 dimensional_split
    # Lin-Rood PPM 2D transport rule (ESS RFC §7.5).
    @test "fv3_lin_rood_advection" in names
    @test "fv3_lin_rood_advection" in Set(r.name for r in fv_rules)
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

@testitem "centered_2nd_deriv_uniform scheme is discoverable and well-formed (esd-8f8)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_deriv_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"d2\"", content)
    # Stencil form with cartesian selectors: three offsets for [1,-2,1]/dx^2.
    @test occursin("\"stencil\"", content)
    @test occursin("\"kind\": \"cartesian\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
    @test occursin("\"offset\": 1", content)
    # No inline replacement — the stencil is lowered by stencil_lowering.jl.
    @test !occursin("\"replacement\"", content)
end

@testitem "laplacian_2nd_uniform_cartesian scheme is discoverable and well-formed (esd-8f8)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "laplacian_2nd_uniform_cartesian", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"laplacian\"", content)
    # Stencil form with arakawa selectors: six entries covering $x and $y axes.
    @test occursin("\"stencil\"", content)
    @test occursin("\"kind\": \"arakawa\"", content)
    @test occursin("\"\$x\"", content)
    @test occursin("\"\$y\"", content)
    @test occursin("\"stagger\": \"cell_center\"", content)
    # All three offsets per axis.
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
    @test occursin("\"offset\": 1", content)
    # No inline replacement — the stencil is lowered by stencil_lowering.jl.
    @test !occursin("\"replacement\"", content)
end

@testitem "nonlinear_laplacian_uniform scheme is discoverable and well-formed (esd-1p7)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "nonlinear_laplacian_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Scheme is identified by a depth-2 grad pattern: Dx(f * Dx(u)).
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Depth-2: the applies_to args contain a nested * with an inner grad.
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"\$f\"", content)
    @test occursin("\"\$u\"", content)
    @test occursin("\"\$x\"", content)
    # Inline replacement — face-interpolated flux difference as arrayop.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # Face interpolation: 0.5*(f[i]+f[i+1]) — must use index and 0.5.
    @test occursin("0.5", content)
    @test occursin("\"op\": \"index\"", content)
    # Flux difference over dx^2.
    @test occursin("\"dx\"", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "spherical_laplacian_uniform scheme is discoverable and well-formed (esd-thf)" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "spherical_laplacian_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Scheme identified by the spherical_laplacian op on the radial dimension.
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"spherical_laplacian\"", content)
    @test occursin("\"\$u\"", content)
    @test occursin("\"\$r\"", content)
    # Inline arrayop replacement — conservative face-flux form.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # Half-point face radii: 0.5*(r[i] + r[i+1]) and 0.5*(r[i-1] + r[i]).
    @test occursin("0.5", content)
    @test occursin("\"op\": \"index\"", content)
    # Coordinate array r and spacing dr appear in denominator.
    @test occursin("\"dr\"", content)
    @test occursin("\"r\"", content)
    # Squared face radii appear as ^ 2.
    @test occursin("\"op\": \"^\"", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_2nd_nonuniform_cartesian scheme is discoverable and well-formed (esd-2t4)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_nonuniform_cartesian", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"d2\"", content)
    # Stencil form with cartesian selectors: three offsets for non-uniform formula.
    @test occursin("\"stencil\"", content)
    @test occursin("\"kind\": \"cartesian\"", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 0", content)
    @test occursin("\"offset\": 1", content)
    # Non-uniform: explicit index nodes for per-cell dx[i] and dx[i+1] coefficients.
    @test occursin("\"index\"", content)
    @test occursin("\"dx\"", content)
    # No inline replacement — the stencil is lowered by stencil_lowering.jl.
    @test !occursin("\"replacement\"", content)
end

@testitem "mixed_deriv_2nd_uniform scheme is discoverable and well-formed (esd-wdv)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "mixed_deriv_2nd_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Depth-2 applies_to: Dx(Dy(u)) — outer grad dim=$x, inner grad dim=$y.
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Nested pattern: inner grad carries dim=$y, outer carries dim=$x.
    @test occursin("\"\$y\"", content)
    @test occursin("\"\$x\"", content)
    @test occursin("\"\$u\"", content)
    # Inline replacement — 4-point stencil as arrayop over 2D output.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # 4-point stencil reads at ±1 offsets in both axes.
    @test occursin("\"op\": \"index\"", content)
    # Denominator: 4*dx*dy.
    @test occursin("\"dx\"", content)
    @test occursin("\"dy\"", content)
    @test occursin("4", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_4th_uniform scheme is discoverable and well-formed (esd-znb)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_4th_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"order\": 4", content)
    # Stencil form with cartesian selectors: four offsets for Fornberg 5-point
    # stencil [1,-8,8,-1]/(12*dx) — offset 0 omitted (zero coefficient).
    @test occursin("\"stencil\"", content)
    @test occursin("\"kind\": \"cartesian\"", content)
    @test occursin("\"offset\": -2", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 1", content)
    @test occursin("\"offset\": 2", content)
    # Fornberg weights: ±1 and ±8 over denominator 12*dx.
    @test occursin("12", content)
    @test occursin("\"dx\"", content)
    # No inline replacement — the stencil is lowered by stencil_lowering.jl.
    @test !occursin("\"replacement\"", content)
    # AST-only: no call op.
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_6th_uniform scheme is discoverable and well-formed (esd-znb)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_6th_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"order\": 6", content)
    # Stencil form with cartesian selectors: six offsets for Fornberg 7-point
    # stencil [-1,9,-45,45,-9,1]/(60*dx) — offset 0 omitted (zero coefficient).
    @test occursin("\"stencil\"", content)
    @test occursin("\"kind\": \"cartesian\"", content)
    @test occursin("\"offset\": -3", content)
    @test occursin("\"offset\": -2", content)
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 1", content)
    @test occursin("\"offset\": 2", content)
    @test occursin("\"offset\": 3", content)
    # Fornberg weights: ±1, ±9, ±45 over denominator 60*dx.
    @test occursin("60", content)
    @test occursin("\"dx\"", content)
    # No inline replacement — the stencil is lowered by stencil_lowering.jl.
    @test !occursin("\"replacement\"", content)
    # AST-only: no call op.
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "mixed_deriv_2nd_nonuniform scheme is discoverable and well-formed (esd-2uy)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "mixed_deriv_2nd_nonuniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Depth-2 applies_to: Dx(Dy(u)) — same pattern as mixed_deriv_2nd_uniform.
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"\$y\"", content)
    @test occursin("\"\$x\"", content)
    @test occursin("\"\$u\"", content)
    # Inline replacement — 4-point stencil as arrayop over 2D output.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"op\": \"index\"", content)
    # Non-uniform denominator: per-cell dx[$x], dx[$x+1], dy[$y], dy[$y+1].
    @test occursin("\"dx\"", content)
    @test occursin("\"dy\"", content)
    # No flat scalar 4*dx*dy denominator — spacing must be indexed per-cell.
    @test !occursin("\"args\": [4, {\"op\": \"*\", \"args\": [\"dx\", \"dy\"]}]", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "nonlinear_laplacian_nonuniform scheme is discoverable and well-formed (esd-2uy)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "nonlinear_laplacian_nonuniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Same depth-2 grad pattern as nonlinear_laplacian_uniform: Dx(f * Dx(u)).
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"\$f\"", content)
    @test occursin("\"\$u\"", content)
    @test occursin("\"\$x\"", content)
    # Inline replacement — conservative flux-difference form as arrayop.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # Face interpolation: 0.5*(f[i]+f[i+1]).
    @test occursin("0.5", content)
    @test occursin("\"op\": \"index\"", content)
    # Non-uniform: per-cell dx[$x] and dx[$x+1] in both gradient and cell-width terms.
    @test occursin("\"dx\"", content)
    # No flat scalar dx² denominator — must be indexed per-cell.
    @test !occursin("\"args\": [\"dx\", \"dx\"]", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "spherical_laplacian_nonuniform scheme is discoverable and well-formed (esd-2uy)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "spherical_laplacian_nonuniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Same applies_to as spherical_laplacian_uniform.
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"spherical_laplacian\"", content)
    @test occursin("\"\$u\"", content)
    @test occursin("\"\$r\"", content)
    # Inline arrayop replacement — conservative face-flux form.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # Half-point face radii: 0.5*(r[i]+r[i+1]) and 0.5*(r[i-1]+r[i]).
    @test occursin("0.5", content)
    @test occursin("\"op\": \"index\"", content)
    # Coordinate array r and per-cell spacing dr appear in denominator.
    @test occursin("\"dr\"", content)
    @test occursin("\"r\"", content)
    # Squared face radii appear as ^ 2.
    @test occursin("\"op\": \"^\"", content)
    # Non-uniform: per-cell dr[$r] and dr[$r+1] instead of flat scalar dr².
    @test !occursin("\"args\": [\"dr\", \"dr\"]", content)
    # No stencil blobs, no call op.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "fv3_lin_rood_advection dimensional_split rule is discoverable and well-formed (esd-f4i)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "fv3_lin_rood_advection", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    # Outer composite rule: dimensional_split discriminator (§7.5).
    @test occursin("\"kind\": \"dimensional_split\"", content)
    # applies_to: 2D advect on cubed-sphere matching transport_2d.json.
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"advect\"", content)
    @test occursin("\"dims\"", content)
    @test occursin("\"xi\"", content)
    @test occursin("\"eta\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cubed_sphere\"", content)
    # Dimensional-split structural fields (§7.5 schema).
    @test occursin("\"axes\"", content)
    @test occursin("\"inner_rule\"", content)
    @test occursin("\"fv3_lin_rood_ppm_1d\"", content)
    @test occursin("\"splitting\"", content)
    @test occursin("\"strang\"", content)
    # Inner rule sibling present in the same file.
    @test occursin("\"fv3_lin_rood_ppm_1d\"", content)
    # No flat stencil or combine on the outer dimensional_split entry (§7.5).
    # (The outer entry has no stencil; stencil may appear inside inner rule docs.)
    # No call op (AGENTS.md constraint).
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "fv3_lin_rood_advection: inner rule is cartesian 1D PPM and outer matches transport_2d shape (esd-f4i)" begin
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(
        repo_root, "discretizations", "finite_volume", "fv3_lin_rood_advection.json"
    )
    @test isfile(rule_path)
    parsed = JSON.parsefile(rule_path)
    discs = parsed["discretizations"]

    # Both sibling entries must be present.
    @test haskey(discs, "fv3_lin_rood_ppm_1d")
    @test haskey(discs, "fv3_lin_rood_advection")

    inner = discs["fv3_lin_rood_ppm_1d"]
    outer = discs["fv3_lin_rood_advection"]

    # Inner rule: 1D PPM, cartesian grid family.
    @test get(inner, "grid_family", nothing) == "cartesian"
    at_inner = inner["applies_to"]
    @test haskey(at_inner, "dim")

    # Outer rule: dimensional_split kind, cubed_sphere, Strang splitting, xi+eta axes.
    @test get(outer, "kind", nothing) == "dimensional_split"
    @test get(outer, "grid_family", nothing) == "cubed_sphere"
    @test get(outer, "splitting", nothing) == "strang"
    axes = outer["axes"]
    @test axes == ["xi", "eta"]
    @test outer["inner_rule"] == "fv3_lin_rood_ppm_1d"

    # applies_to matches the transport_2d.json op/args shape: advect($q,$U) over dims.
    at = outer["applies_to"]
    @test get(at, "op", nothing) == "advect"
    @test haskey(at, "args")
    @test haskey(at, "dims")
end

@testitem "fv3_lin_rood_advection: sweep consumes declaration — fields map to transport_2d_linrood! (esd-f4i)" begin
    using EarthSciDiscretizations
    using JSON

    # Confirm that the dimensional_split declaration's fields correspond exactly
    # to what transport_2d_linrood! implements (ESD sweep loop owns execution).
    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(
        repo_root, "discretizations", "finite_volume", "fv3_lin_rood_advection.json"
    )
    outer = JSON.parsefile(rule_path)["discretizations"]["fv3_lin_rood_advection"]

    # axes: the sweep iterates over xi then eta (cubed-sphere in-panel directions).
    @test outer["axes"] == ["xi", "eta"]

    # splitting: "strang" — Lin-Rood uses symmetric half-step advective
    # correction + full flux step (transport_2d_linrood! step 1-3), achieving
    # O(dt^2) accuracy like Strang splitting.
    @test outer["splitting"] == "strang"

    # inner_rule: the 1D operator applied per axis is PPM (flux_1d_ppm! in
    # the sweep; fv3_lin_rood_ppm_1d in the declaration).
    @test outer["inner_rule"] == "fv3_lin_rood_ppm_1d"

    # Functional contract: constant field with zero velocity has zero tendency
    # everywhere (transport_2d_linrood! exact null with no flow; matching the
    # test in test_transport_2d.jl::Lin-Rood 2D transport of constant field).
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)
    q = fill(3.14, 6, Nc, Nc)
    vel_xi  = fill(0.0, 6, Nc + 1, Nc)
    vel_eta = fill(0.0, 6, Nc, Nc + 1)
    dt = 0.01
    tendency = zeros(6, Nc, Nc)
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)
    @test all(abs.(tendency) .< 1.0e-12)
end
