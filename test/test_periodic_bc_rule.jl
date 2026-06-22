using Test
using TestItems

# ===========================================================================
# periodic_bc rule (esd-7mj): the periodic wrap as a declarative makearray-region
# BC rule in the LANDED neumann/robin encoding (fn/dim `bc` wrapper + ess-hjg
# makearray ghost splice), authored fresh to SUPERSEDE the stale esd-agh
# symbolic-`Nx`/`mod` index-rewrite rule.
#
# The ghost reads the wrapped OPPOSITE-end cell via a single side-generic rule:
# the 0-based offset L = N-1 (N = side axis size, bound by the EXISTING
# `bind_side_dim_size` guard) is re-indexed by the makearray lowering
# (`_reindex_ghost`: min -> 1+L, max -> N-L) to land on
#   xmin: 1+(N-1) = N  -> u[N]  (the i-1 ghost of the first cell reads the last)
#   xmax: N-(N-1) = 1  -> u[1]  (the i+1 ghost of the last cell reads the first)
# closing the axis into a torus with NO symbolic Nx, NO mod, NO new primitive.
# ===========================================================================

# ---------------------------------------------------------------------------
# Rule UNIT: the rule fires side-generically and emits index(u, N-1) for each
# axis side, with N concrete from the grid (NO mod / NO Nx in the output).
# ---------------------------------------------------------------------------
@testitem "periodic_bc rule: side-generic opposite-end ghost index(u, N-1) (esd-7mj)" tags = [:bc, :periodic] begin
    include(joinpath(@__DIR__, "walk_esd_tests.jl"))
    using .WalkESDTests
    using EarthSciDiscretizations
    import EarthSciSerialization
    import JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_difference", "periodic_bc.json")
    rules_obj = JSON.parse(read(rule_path, String))["rules"]
    rules = EarthSciSerialization.parse_rules(rules_obj)

    function ghost(side, dim_sizes, spatial)
        ctx = WalkESDTests._build_rule_context(
            Dict{String, Any}(
                "grids" => Dict{String, Any}(
                    "g1" => Dict{String, Any}(
                        "spatial_dims" => spatial, "dim_sizes" => dim_sizes
                    )
                ),
                "variables" => Dict{String, Any}("u" => Dict{String, Any}("grid" => "g1")),
            )
        )
        bc = Dict{String, Any}("op" => "bc", "fn" => "periodic", "dim" => side, "args" => ["u"])
        out = EarthSciSerialization.rewrite(WalkESDTests._expr_from_json(bc), rules, ctx)
        return EarthSciSerialization.canonical_json(out)
    end

    # Axis 1 (x), N=4: both sides author the SAME 0-based offset L = N-1 = 4-1.
    # The makearray lowering re-indexes it to the opposite end per side.
    xmin = ghost("xmin", Dict{String, Any}("x" => 4), ["x"])
    xmax = ghost("xmax", Dict{String, Any}("x" => 4), ["x"])
    @test xmin == "{\"args\":[\"u\",{\"args\":[4,1],\"op\":\"-\"}],\"op\":\"index\"}"
    @test xmax == xmin   # one side-generic rule -> identical 0-based ghost

    # Axis 2 (y), Ny=6: the side-generic `bind_side_dim_size` resolves axis 2's
    # size, proving the rule is axis-position-generic (covers y / lon, not just x).
    ymin = ghost("ymin", Dict{String, Any}("x" => 4, "y" => 6), ["x", "y"])
    ymax = ghost("ymax", Dict{String, Any}("x" => 4, "y" => 6), ["x", "y"])
    @test ymin == "{\"args\":[\"u\",{\"args\":[6,1],\"op\":\"-\"}],\"op\":\"index\"}"
    @test ymax == ymin

    # DECLARATIVE-OR-FAIL contract: the emitted ghost is concrete (N folded in by
    # the bind_side_dim_size binding) — never the retired symbolic Nx / mod form.
    for g in (xmin, xmax, ymin, ymax)
        @test !occursin("mod", g)
        @test !occursin("Nx", g)
    end
end

# ---------------------------------------------------------------------------
# Per-axis SMOKE (1D): build_ode_problem with model-level `kind:"periodic"` BCs
# on a NON-periodic grid -> the makearray lowering splices the opposite-end
# ghost, so the centered-FD Laplacian wraps (du[1] reads u[N], du[N] reads u[1]).
# ---------------------------------------------------------------------------
@testitem "periodic_bc integration: 1D centered-FD Laplacian wraps via makearray ghost (esd-7mj)" tags = [:bc, :periodic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "bc_periodic_1d.esm")
    gdd_path = joinpath(repo_root, "discretizations", "gdd", "cartesian_1d_nonperiodic_n8.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
    N = 8; dx = 1.0 / N
    @test length(prob.u0) == N

    u = [sin(2π * (i - 0.5) * dx) for i in 1:N]
    u0 = copy(prob.u0)
    for i in 1:N
        u0[var_map["u[$i]"]] = u[i]
    end
    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    # Periodic centered second-derivative: du[i] = (u[i-1] - 2u[i] + u[i+1])/dx²
    # with i-1=0 -> u[N] (xmin ghost) and i+1=N+1 -> u[1] (xmax ghost).
    for i in 1:N
        im1 = mod1(i - 1, N); ip1 = mod1(i + 1, N)
        expected = (u[im1] - 2 * u[i] + u[ip1]) / dx^2
        @test du[var_map["u[$i]"]] ≈ expected atol = 1.0e-10
    end
    # Guard the wrap explicitly: boundary cells MUST see the opposite end.
    @test du[var_map["u[1]"]] ≈ (u[N] - 2 * u[1] + u[2]) / dx^2 atol = 1.0e-10
    @test du[var_map["u[$N]"]] ≈ (u[N - 1] - 2 * u[N] + u[1]) / dx^2 atol = 1.0e-10
end

# ---------------------------------------------------------------------------
# Per-axis SMOKE (2D + CORNERS): both x and y periodic. The 5-point Laplacian
# wraps on both axes, corners included. This exercises the makearray ghost path
# end to end: the `laplacian_2nd_uniform_cartesian` rule authors additive offsets
# `index(u, i + (-1), j)`; ESS canonicalize reorders the commutative `+` to
# constant-first `[-1, "i"]`. ESS `_scan_stencil_reach!` now detects that form via
# its two-arm scan (bead ess-wg0). Before that fix the per-axis reach scanned to
# 0, `_apply_makearray_bcs!` emitted NO boundary regions, and every boundary read
# fell back to the zero-ghost convention (the 1-D `d2` rule escaped only because
# it authors one neighbour as a non-commutative SUBTRACTION). With ess-wg0 landed
# the boundary/corner reads splice the periodic ghost, so the wrap is byte-exact
# against the grid-periodic ground truth — see discretizations/finite_difference/
# periodic_bc/MAKEARRAY_2D_REACH_GAP.md.
# ---------------------------------------------------------------------------
@testitem "periodic_bc integration: 2D 5-point Laplacian wraps both axes + corners (esd-7mj)" tags = [:bc, :periodic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "bc_periodic_2d.esm")
    gdd_path = joinpath(repo_root, "discretizations", "gdd", "cartesian_2d_nonperiodic_n8.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
    N = 8; dx = 1.0 / N; dy = 1.0 / N
    @test length(prob.u0) == N * N

    # Asymmetric IC so the periodic wrap differs from the zero-ghost fallback at
    # EVERY boundary cell (a symmetric sin·cos can make a corner coincide).
    U = zeros(N, N)
    u0 = copy(prob.u0)
    for i in 1:N, j in 1:N
        v = sin(2π * (i - 0.5) * dx) * cos(2π * (j - 0.5) * dy) + 0.4 * cos(2π * (i - 0.5) * dx)
        U[i, j] = v
        u0[var_map["u[$i,$j]"]] = v
    end
    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    analytic(i, j) = begin
        il = mod1(i - 1, N); ir = mod1(i + 1, N)
        jb = mod1(j - 1, N); ja = mod1(j + 1, N)
        (U[il, j] - 2 * U[i, j] + U[ir, j]) / dx^2 +
            (U[i, jb] - 2 * U[i, j] + U[i, ja]) / dy^2
    end

    # Interior cells read no out-of-range neighbour, so they match TODAY
    # (region 0 carries the raw stencil) — independent of the 2-D reach gap.
    for i in 2:(N - 1), j in 2:(N - 1)
        @test du[var_map["u[$i,$j]"]] ≈ analytic(i, j) atol = 1.0e-10
    end

    # Boundary + corner wrap: enabled by the ESS _scan_stencil_reach! two-arm scan
    # (ess-wg0). The periodic rule fires (all four BCs receive index(u, N-1)) and
    # the makearray box-shrink now keeps the boundary regions, so the ghost splices.
    boundary = [(i, j) for i in 1:N for j in 1:N if i == 1 || i == N || j == 1 || j == N]
    max_boundary_err = maximum(abs(du[var_map["u[$i,$j]"]] - analytic(i, j)) for (i, j) in boundary)
    @test max_boundary_err < 1.0e-10
    # Corner (1,1) must wrap on BOTH axes (x-neighbour u[N,1], y-neighbour u[1,N]).
    @test du[var_map["u[1,1]"]] ≈ analytic(1, 1) atol = 1.0e-10
end
