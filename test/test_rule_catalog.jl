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
    @test occursin("\"offset\": -1", content)
    @test occursin("\"offset\": 1", content)
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
        "periodic_bc",
        "upwind_1st",
    )
        @test seeded in names
    end
    fd_rules = filter(r -> r.family == :finite_difference, rules)
    fd_names = Set(r.name for r in fd_rules)
    @test "centered_2nd_uniform" in fd_names
    @test "centered_2nd_uniform_vertical" in fd_names
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
end
