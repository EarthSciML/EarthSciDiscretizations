using Test
using TestItems

# BC/IC golden tests (esd-7i3).
# These pin byte-identity for the bc-ic-unification refactor campaign.
# Supported cases (2-D Dirichlet, IC uniform/nonuniform) run live.
# Unsupported cases (nonzero-Neumann, Robin): @test_throws pins current behavior;
# analytic expected values are computed in the same testitem for future reference.

# ---------------------------------------------------------------------------
# IC goldens: pin _eval_expression_ics output
# ---------------------------------------------------------------------------

@testitem "bc_ic goldens: IC sin(2π·x) on uniform 1D N=8 (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "ic_uniform_1d.esm")
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "cartesian_1d_nonperiodic_n8.gdd.json"
    )

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N = 8; dx = 0.125
    @test length(prob.u0) == N
    for i in 1:N
        x = (i - 0.5) * dx
        @test prob.u0[var_map["u[$i]"]] ≈ sin(2π * x)  atol = 1e-15
    end
end

@testitem "bc_ic goldens: IC sin(2π·z) on nonuniform vertical N=16 (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "ic_nonuniform_1d.esm")
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "vertical_1d_nonuniform_n16.gdd.json"
    )

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    # Reproduce cell centres as _eval_expression_ics does (cumsum(widths) - widths/2).
    levels = [
        0.0, 0.07468119, 0.14750791, 0.216908, 0.28183099, 0.341908,
        0.39750791, 0.44968119, 0.5, 0.55031881, 0.60249209, 0.658092,
        0.71816901, 0.783092, 0.85249209, 0.92531881, 1.0,
    ]
    widths = diff(levels)
    N = length(widths)
    cumw = cumsum(widths)
    centres = [cumw[k] - widths[k] / 2 for k in 1:N]

    @test length(prob.u0) >= N
    for i in 1:N
        @test prob.u0[var_map["u[$i]"]] ≈ sin(2π * centres[i])  atol = 1e-14
    end
end

# ---------------------------------------------------------------------------
# 2-D Dirichlet non-periodic: pin du byte-identity (all four sides + corners)
# ---------------------------------------------------------------------------

@testitem "bc_ic goldens: 2D Dirichlet-0 nonperiodic Laplacian du (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(
        repo_root, "test", "fixtures", "bc_ic", "bc_dirichlet_2d_nonperiodic.esm"
    )
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "cartesian_2d_nonperiodic_n8.gdd.json"
    )

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N = 8; dx = 0.125; dy = 0.125
    @test length(prob.u0) == N * N

    # IC: sin(π·x)·sin(π·y) at cell centres (ESM ic expression)
    u0 = copy(prob.u0)
    u = zeros(N, N)
    for i in 1:N, j in 1:N
        x = (i - 0.5) * dx; y = (j - 0.5) * dy
        v = sin(π * x) * sin(π * y)
        u[i, j] = v
        u0[var_map["u[$i,$j]"]] = v
    end

    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    # Expected: 5-point Laplacian with Dirichlet-0 ghost (ghost = bc value = 0)
    for i in 1:N, j in 1:N
        u_left  = i > 1 ? u[i-1, j] : 0.0
        u_right = i < N ? u[i+1, j] : 0.0
        u_below = j > 1 ? u[i, j-1] : 0.0
        u_above = j < N ? u[i, j+1] : 0.0
        expected = (u_left - 2*u[i,j] + u_right) / dx^2 +
                   (u_below - 2*u[i,j] + u_above) / dy^2
        @test du[var_map["u[$i,$j]"]] ≈ expected  atol = 1e-10
    end
end

# ---------------------------------------------------------------------------
# 2-D periodic-x + Dirichlet-y: pin du byte-identity
# ---------------------------------------------------------------------------

@testitem "bc_ic goldens: 2D periodic-x + Dirichlet-y Laplacian du (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(
        repo_root, "test", "fixtures", "bc_ic", "bc_dirichlet_2d_periodic_x.esm"
    )
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd",
        "cartesian_2d_periodic_x_dirichlet_y_n8.gdd.json",
    )

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N = 8; dx = 0.125; dy = 0.125
    @test length(prob.u0) == N * N

    # IC: sin(2π·x)·sin(π·y) at cell centres (ESM ic expression)
    u0 = copy(prob.u0)
    u = zeros(N, N)
    for i in 1:N, j in 1:N
        x = (i - 0.5) * dx; y = (j - 0.5) * dy
        v = sin(2π * x) * sin(π * y)
        u[i, j] = v
        u0[var_map["u[$i,$j]"]] = v
    end

    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    # Expected: x is periodic (mod-wrap), y has Dirichlet-0 ghosts
    for i in 1:N, j in 1:N
        u_left  = u[mod1(i - 1, N), j]       # periodic x
        u_right = u[mod1(i + 1, N), j]       # periodic x
        u_below = j > 1 ? u[i, j-1] : 0.0   # Dirichlet ymin = 0
        u_above = j < N ? u[i, j+1] : 0.0   # Dirichlet ymax = 0
        expected = (u_left - 2*u[i,j] + u_right) / dx^2 +
                   (u_below - 2*u[i,j] + u_above) / dy^2
        @test du[var_map["u[$i,$j]"]] ≈ expected  atol = 1e-10
    end
end

# ---------------------------------------------------------------------------
# Nonzero-Neumann 1D: IC-arrayop now succeeds; du uses ESS default ghost = 0
# ---------------------------------------------------------------------------

@testitem "bc_ic goldens: nonzero-Neumann 1D golden u0 + default-ghost du (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "bc_neumann_1d_nonzero.esm")
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "cartesian_1d_nonperiodic_n8.gdd.json"
    )

    # _prepare_expression_ics! converts type:expression IC to initialization_equations before
    # calling ESS discretize, so build_ode_problem now succeeds without throwing.
    # The neumann_bc.json rule is not wired into this GDD's discretization set; ESS defaults
    # to ghost = 0 (Dirichlet-0) for out-of-bounds indices instead of the Neumann prescription.
    # Analytic golden (ghost formula from neumann_bc.json, for future GDD integration):
    #   u_ghost = u[boundary_cell] + dx * value
    #   xmin BC: value=1.0 → ghost_xmin = u[1] + dx  (not yet applied)
    #   xmax BC: value=0.0 → ghost_xmax = u[N]        (not yet applied)
    N = 8; dx = 0.125
    u = [sin(π * (i - 0.5) * dx) for i in 1:N]
    # ESS uses ghost = 0 (Dirichlet-0 default) — Neumann rule not in GDD discretizations.
    ghost_xmin = 0.0
    ghost_xmax = 0.0
    expected = zeros(N)
    expected[1] = (ghost_xmin - 2*u[1] + u[2]) / dx^2
    expected[N] = (u[N-1] - 2*u[N] + ghost_xmax) / dx^2
    for i in 2:N-1
        expected[i] = (u[i-1] - 2*u[i] + u[i+1]) / dx^2
    end

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
    @test length(prob.u0) == N
    for i in 1:N
        @test prob.u0[var_map["u[$i]"]] ≈ u[i]  atol = 1e-14
    end
    du = similar(prob.u0)
    prob.f(du, prob.u0, prob.p, 0.0)
    for i in 1:N
        @test du[var_map["u[$i]"]] ≈ expected[i]  atol = 1e-10
    end
end

# ---------------------------------------------------------------------------
# Robin 1D: IC-arrayop now succeeds; du uses ESS default ghost = 0
# ---------------------------------------------------------------------------

@testitem "bc_ic goldens: Robin 1D golden u0 + default-ghost du (esd-7i3)" tags=[:bc_ic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "bc_ic", "bc_robin_1d.esm")
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "cartesian_1d_nonperiodic_n8.gdd.json"
    )

    # _prepare_expression_ics! converts type:expression IC to initialization_equations before
    # calling ESS discretize, so build_ode_problem now succeeds without throwing.
    # The robin_bc.json rule is not wired into this GDD's discretization set; ESS defaults
    # to ghost = 0 (Dirichlet-0) for out-of-bounds indices instead of the Robin prescription.
    # Analytic golden (ghost formula from robin_bc.json, for future GDD integration):
    #   u_ghost = (2*dx*gamma + (2*beta - alpha*dx)*u[1]) / (alpha*dx + 2*beta)
    #   alpha=1, beta=1, gamma=0.5, dx=0.125: ghost_xmin = (0.125 + 1.875*u[1]) / 2.125
    #   xmax: no BC declared → ghost_xmax = 0
    N = 8; dx = 0.125
    u = [sin(π * (i - 0.5) * dx) for i in 1:N]
    # ESS uses ghost = 0 (Dirichlet-0 default) — Robin rule not in GDD discretizations.
    ghost_xmin = 0.0
    ghost_xmax = 0.0
    expected = zeros(N)
    expected[1] = (ghost_xmin - 2*u[1] + u[2]) / dx^2
    expected[N] = (u[N-1] - 2*u[N] + ghost_xmax) / dx^2
    for i in 2:N-1
        expected[i] = (u[i-1] - 2*u[i] + u[i+1]) / dx^2
    end

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
    @test length(prob.u0) == N
    for i in 1:N
        @test prob.u0[var_map["u[$i]"]] ≈ u[i]  atol = 1e-14
    end
    du = similar(prob.u0)
    prob.f(du, prob.u0, prob.p, 0.0)
    for i in 1:N
        @test du[var_map["u[$i]"]] ≈ expected[i]  atol = 1e-10
    end
end
