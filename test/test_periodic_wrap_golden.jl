using Test
using TestItems

# ===========================================================================
# Cross-rig periodic regression GOLDEN (esd-2my).
#
# ess-8ne (ESS d9049910) retired the imperative `_apply_periodic_folding!` fold:
# periodic discretization now flows ENTIRELY through the declarative
# makearray-region ghost path (the same path ess-hjg built for declarative BCs).
# The 2-D boundary/corner ghost only splices once `_scan_stencil_reach!` detects
# the constant-first stencil offset the `laplacian_2nd_uniform_cartesian` rule
# authors after ESS canonicalize reorders the commutative `+` (bead ess-wg0).
#
# This golden exists because that 2-D reach bug went unnoticed for a full
# refinery cycle: NO ESD end-to-end test exercised a non-trivial periodic ghost
# coming from the GRID itself (grid-level `boundary_conditions: periodic`, as
# opposed to the model-level `kind:"periodic"` opt-in already pinned in
# test_periodic_bc_rule.jl). It closes that coverage gap and makes the
# cross-rig verification a CI-gated commit:
#
#   * the 5-point Laplacian built from a grid-level-periodic GDD must wrap on
#     BOTH axes and at the corners, reproducing the analytic periodic operator
#     to machine precision, and
#   * the wrapped corner value is byte-pinned far from the zero-ghost fallback,
#     so a silent reach/fold regression (du[1,1] would jump from ≈-30 to ≈-231)
#     fails this test loudly rather than passing on the interior alone.
# ===========================================================================

@testitem "periodic golden: 2D grid-level periodic Laplacian wraps both axes + corners (esd-2my)" tags = [:bc, :periodic] begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    # Grid-level periodicity: the GDD declares `boundary_conditions: [{type:
    # periodic, dimensions: [x, y]}]` — the periodicity lives in the grid, not
    # in a model-level BC opt-in. diffusion_2d_pde.esm carries D(u,t)=laplacian(u).
    esm_path = joinpath(
        repo_root, "discretizations", "finite_difference",
        "laplacian_2nd_uniform_cartesian", "fixtures", "integration",
        "diffusion_2d_pde.esm",
    )
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd", "cartesian_2d_periodic_n16.gdd.json"
    )

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
    N = 16; dx = 1.0 / N; dy = 1.0 / N
    @test length(prob.u0) == N * N

    # Asymmetric IC so the periodic wrap differs from the zero-ghost fallback at
    # EVERY boundary cell (a symmetric sin·cos can make a corner coincide). We
    # overwrite ALL N² state entries, so the .esm's own sin·sin IC drops out —
    # this isolates the discretized operator, not the IC.
    U = zeros(N, N)
    u0 = copy(prob.u0)
    for i in 1:N, j in 1:N
        v = sin(2π * (i - 0.5) * dx) * cos(2π * (j - 0.5) * dy) + 0.4 * cos(2π * (i - 0.5) * dx)
        U[i, j] = v
        u0[var_map["u[$i,$j]"]] = v
    end
    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    # Analytic periodic 5-point Laplacian (mod1 closes each axis into a torus).
    analytic(i, j) = begin
        il = mod1(i - 1, N); ir = mod1(i + 1, N)
        jb = mod1(j - 1, N); ja = mod1(j + 1, N)
        (U[il, j] - 2 * U[i, j] + U[ir, j]) / dx^2 +
            (U[i, jb] - 2 * U[i, j] + U[i, ja]) / dy^2
    end

    # Whole field — interior AND every boundary/corner cell — reproduces the
    # analytic PERIODIC operator. The boundary arm is the cross-rig regression
    # check: it only holds because the grid-level periodic ghost splices through
    # the ess-wg0 / ess-8ne makearray-region path.
    max_err = maximum(abs(du[var_map["u[$i,$j]"]] - analytic(i, j)) for i in 1:N for j in 1:N)
    @test max_err < 1.0e-9

    # Corner (1,1) must wrap on BOTH axes: x-neighbour u[N,1], y-neighbour u[1,N].
    @test du[var_map["u[1,1]"]] ≈ analytic(1, 1) atol = 1.0e-9
    @test du[var_map["u[$N,$N]"]] ≈ analytic(N, N) atol = 1.0e-9

    # Byte-pinned golden + non-triviality guard. The wrapped corner value is
    # ≈-30.20; the NON-wrapping (zero-ghost) value would be ≈-231.07. Freezing
    # the wrapped value and asserting it is nowhere near the zero-ghost fallback
    # makes a silent un-wrap regression impossible to pass.
    golden_du_11 = -30.20450580068021
    @test du[var_map["u[1,1]"]] ≈ golden_du_11 rtol = 1.0e-10
    zero_ghost_du_11 = (
        (0.0 - 2 * U[1, 1] + U[2, 1]) / dx^2 +
            (0.0 - 2 * U[1, 1] + U[1, 2]) / dy^2
    )
    @test abs(du[var_map["u[1,1]"]] - zero_ghost_du_11) > 100.0
end
