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
    @test occursin("\"op\": \"grad\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields retired. The replacement encodes (u[$x] - u[$x-1]) / dx.
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
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
    @test occursin("\"op\": \"grad\"", content)
    # Migrated to arrayop replacement form (esd-7h2): no longer stencil.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
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
    # `kind:"dirichlet"` and `side:"$side"` (generalized from xmin — esd-klj),
    # `replacement` giving the ghost-cell value 2*u_bc - u[0], and `produces`
    # declaring the ghost_var per §9.4. No $h needed: ghost is spacing-independent.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"dirichlet\"", content)
    @test occursin("\"side\": \"\$side\"", content)
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
    @test occursin("__\$side", content)
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
    # `kind:"neumann"` and `side:"$side"` (generalized from xmin — esd-klj),
    # `replacement` giving the ghost-cell value u[0] + $h*value (where value is
    # du/dn, the outward normal derivative, and $h is the grid spacing for the
    # side's axis, bound via ESS bind_side_spacing), and `produces` per §9.4.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"neumann\"", content)
    @test occursin("\"side\": \"\$side\"", content)
    @test occursin("\"op\": \"bc\"", content)
    # Ghost-cell formula: u_ghost = u[0] + $h*value.
    @test occursin("\"op\": \"+\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"\$h\"", content)
    # §9.4 ghost_var produces declaration.
    @test occursin("\"produces\"", content)
    @test occursin("\"ghost_var\"", content)
    # Ghost variable named per §9.4 scheme__logical__side convention.
    @test occursin("neumann_bc__", content)
    @test occursin("__\$side", content)
end

@testitem "robin_bc rule is discoverable and well-formed (esd-m9v, generalized esd-klj)" begin
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
    # `kind:"robin"` and `side:"$side"` (generalized from xmin — esd-klj),
    # `robin_alpha/beta/gamma` coefficients for αu + β∂u/∂n = γ,
    # `replacement` giving the ghost-cell value
    # u_ghost = (2·$h·γ + (2·β - α·$h)·u[0]) / (α·$h + 2·β),
    # where $h is the grid spacing for the side's axis (bound via ESS
    # bind_side_spacing), and `produces` declaring the ghost_var per §9.4.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"kind\": \"robin\"", content)
    @test occursin("\"side\": \"\$side\"", content)
    @test occursin("\"op\": \"bc\"", content)
    # Coefficients: robin_alpha (α), robin_beta (β), robin_gamma (γ).
    @test occursin("\"robin_alpha\"", content)
    @test occursin("\"robin_beta\"", content)
    @test occursin("\"robin_gamma\"", content)
    # Ghost-cell formula uses division, addition, subtraction, multiplication,
    # index into u[0], and the side-axis spacing $h.
    @test occursin("\"op\": \"/\"", content)
    @test occursin("\"op\": \"+\"", content)
    @test occursin("\"op\": \"-\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"\$h\"", content)
    # §9.4 ghost_var produces declaration.
    @test occursin("\"produces\"", content)
    @test occursin("\"ghost_var\"", content)
    # Ghost variable named per §9.4 scheme__logical__side convention.
    @test occursin("robin_bc__", content)
    @test occursin("__\$side", content)
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
    # esd-t4h: migrated from stencil form to authored arrayop replacement.
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Replacement uses the pattern variable from applies_to.dim and scalar spacing h.
    @test occursin("\"\$k\"", content)
    @test occursin("\"h\"", content)
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
    # esd-t4h: migrated from stencil form to authored arrayop replacement.
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # Replacement indexes $u along "lat" and "lon" axes with ±1 offsets.
    @test occursin("\"lon\"", content)
    @test occursin("\"lat\"", content)
    # Coefficient symbols: angular spacings dlon/dlat, sphere radius R, and
    # the latitude metric cos_lat (lon-axis only).
    @test occursin("\"dlon\"", content)
    @test occursin("\"dlat\"", content)
    @test occursin("\"R\"", content)
    @test occursin("\"cos_lat\"", content)
end

@testitem "nn_diffusion_duo scheme is discoverable and well-formed" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "nn_diffusion_duo", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"unstructured\"", content)
    # esd-agh: rule migrated from stencil form to authored reduce-arrayop replacement.
    @test occursin("\"replacement\"", content)
    @test !occursin("\"stencil\"", content)
    # Operator class: scalar Laplacian acting at cell centers.
    @test occursin("\"op\": \"laplacian\"", content)
    @test occursin("\"emits_location\": \"cell_center\"", content)
    # Replacement: reduce-arrayop with constant 1-based k∈[1,3] (DUO always has 3 neighbors).
    @test occursin("\"reduce\": \"+\"", content)
    @test occursin("\"ranges\"", content)
    # Coefficient/connectivity symbols: dc_edge, dv_edge, area, edges_on_face, cell_neighbors.
    @test occursin("\"dc_edge\"", content)
    @test occursin("\"dv_edge\"", content)
    @test occursin("\"area\"", content)
    @test occursin("\"cell_neighbors\"", content)
    @test occursin("\"edges_on_face\"", content)
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
    # esd-agh: rule migrated from stencil form to authored reduce-arrayop replacement.
    @test occursin("\"replacement\"", content)
    @test !occursin("\"stencil\"", content)
    # Operator class: scalar Laplacian acting at cell centers.
    @test occursin("\"op\": \"laplacian\"", content)
    @test occursin("\"emits_location\": \"cell_center\"", content)
    # Replacement: reduce-arrayop with variable 1-based range k∈[1, n_edges_on_cell[i]].
    @test occursin("\"reduce\": \"+\"", content)
    @test occursin("\"ranges\"", content)
    @test occursin("\"arrayop\"", content)
    # Coefficient symbols (MPAS snake_case): area_cell, dc_edge, dv_edge,
    # edges_on_cell, cells_on_cell, n_edges_on_cell.
    @test occursin("\"dv_edge\"", content)
    @test occursin("\"dc_edge\"", content)
    @test occursin("\"area_cell\"", content)
    @test occursin("\"edges_on_cell\"", content)
    @test occursin("\"cells_on_cell\"", content)
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
            "nn_diffusion_duo",
            "periodic_bc",
            "upwind_1st",
            "upwind_1st_nonuniform",
            "dirichlet_bc",
            "neumann_bc",
            "robin_bc",
            "mixed_deriv_2nd_nonuniform",
            "nonlinear_laplacian_nonuniform",
            "spherical_laplacian_nonuniform",
            "expression_ic",
        )
        @test seeded in names
    end
    fd_rules = filter(r -> r.family == :finite_difference, rules)
    fd_names = Set(r.name for r in fd_rules)
    @test "centered_2nd_uniform" in fd_names
    @test "centered_2nd_uniform_vertical" in fd_names
    @test "centered_2nd_uniform_latlon" in fd_names
    @test "nn_diffusion_mpas" in fd_names
    @test "nn_diffusion_duo" in fd_names
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
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (u[i-1] + u[i+1] - 2*u[i]) / (dx*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
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
    # Canonical bare replacement (esd-eg5): stencil/arrayop retired.
    # Six terms: u[i±1,j]/dx^2 and u[i,j±1]/dy^2 with canonical i/j indices.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"arrayop\"", content)
    @test !occursin("\"output_idx\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"i\"", content)
    @test occursin("\"j\"", content)
    @test occursin("\"dx\"", content)
    @test occursin("\"dy\"", content)
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
    # Migrated to arrayop replacement form (esd-7h2): no longer stencil.
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    # Non-uniform: explicit index nodes for per-cell dx[i] and dx[i+1] coefficients.
    @test occursin("\"index\"", content)
    @test occursin("\"dx\"", content)
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
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (8*(u[i+1]-u[i-1]) - (u[i+2]-u[i-2])) / (12*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("12", content)
    @test occursin("\"dx\"", content)
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
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (45*(u[i+1]-u[i-1]) - 9*(u[i+2]-u[i-2]) + (u[i+3]-u[i-3])) / (60*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("60", content)
    @test occursin("\"dx\"", content)
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

@testitem "centered_2nd_nonuniform_vertical scheme is discoverable and well-formed (esd-9qs)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_2nd_nonuniform_vertical", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"vertical\"", content)
    @test occursin("\"op\": \"grad\"", content)
    # esd-t4h: migrated from stencil form to authored arrayop replacement.
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    # Per-cell dz[k] binding via index op — the S2 per-cell binding contract (SELECTOR_KINDS #3).
    @test occursin("\"op\": \"index\"", content)
    @test occursin("\"dz\"", content)
    @test occursin("\"\$k\"", content)
    # No flat scalar spacing — must always be indexed.
    @test !occursin("\"args\": [-1, \"dz\"]", content)
    @test !occursin("\"args\": [1, \"dz\"]", content)
    # No call op.
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_8th_uniform scheme is discoverable and well-formed (esd-ec4)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_8th_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (672*(u[i+1]-u[i-1]) - 168*(u[i+2]-u[i-2]) + 32*(u[i+3]-u[i-3]) - 3*(u[i+4]-u[i-4])) / (840*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("840", content)
    @test occursin("\"dx\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_4th_deriv_uniform scheme is discoverable and well-formed (esd-ec4)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_4th_deriv_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"d2\"", content)
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (16*(u[i-1]+u[i+1]) - (u[i-2]+u[i+2]) - 30*u[i]) / (12*dx*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("12", content)
    @test occursin("\"dx\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_6th_deriv_uniform scheme is discoverable and well-formed (esd-ec4)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_6th_deriv_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"d2\"", content)
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (270*(u[i-1]+u[i+1]) - 27*(u[i-2]+u[i+2]) + 2*(u[i-3]+u[i+3]) - 490*u[i]) / (180*dx*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("180", content)
    @test occursin("\"dx\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "centered_8th_deriv_uniform scheme is discoverable and well-formed (esd-ec4)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "centered_8th_deriv_uniform", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"op\": \"d2\"", content)
    @test occursin("\"accuracy\"", content)
    # esd-3d7: migrated to arrayop replacement form — stencil/selector/offset
    # fields and fornberg_gen.py header retired.
    # The replacement encodes (8064*(u[i-1]+u[i+1]) - 1008*(u[i-2]+u[i+2]) + 128*(u[i-3]+u[i+3]) - 9*(u[i-4]+u[i+4]) - 14350*u[i]) / (5040*dx*dx).
    @test !occursin("\"stencil\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"op\": \"index\"", content)
    @test occursin("5040", content)
    @test occursin("\"dx\"", content)
    @test !occursin("\"op\": \"call\"", content)
end

@testitem "expression_ic rule is discoverable and well-formed (esd-klj)" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "expression_ic", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :ic
    @test isfile(rule.path)

    content = read(rule.path, String)
    # expression_ic is a rewrite rule that transforms an `ic` synthetic node
    # into an arrayop over the grid (esd-klj). The engine substitutes spatial-dim
    # symbols with cell-centre coordinates (x→coord_x($i), etc.). ghost_width:0
    # declares that IC evaluation consumes no ghost cells.
    @test occursin("\"pattern\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"op\": \"ic\"", content)
    @test occursin("\"op\": \"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"\$expr\"", content)
    @test occursin("\"ghost_width\"", content)
    @test occursin("0", content)
end
