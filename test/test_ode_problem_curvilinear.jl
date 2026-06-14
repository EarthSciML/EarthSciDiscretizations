using Test
using TestItems
import SciMLBase

# ---------------------------------------------------------------------------
# Path B: build_ode_problem routes curvilinear GDDs through the ESS
# PDESystem → discretize(sys, grid) pipeline (esd-ncg).
# ---------------------------------------------------------------------------

@testitem "build_ode_problem Path-B: LatLon returns solvable ODEProblem" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "latlon_diffusion.esm")
    gdd_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "latlon_diffusion.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem
    @test var_map isa Dict{String, Int}

    # GDD specifies nlon=6, nlat=6 → n_cells = 36
    @test length(prob.u0) == 36

    # Smoke: f!(du, u0, p, t) must execute without error
    du = similar(prob.u0)
    @test_nowarn prob.f(du, prob.u0, prob.p, 0.0)
    @test all(isfinite, du)
end

@testitem "build_ode_problem Path-B: CubedSphere smoke test" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "cubed_sphere_smoke.esm")
    gdd_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "cubed_sphere_smoke.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem
    @test var_map isa Dict{String, Int}

    # GDD: Nc=4 → n_cells = 6 * 4 * 4 = 96
    @test length(prob.u0) == 96

    # Smoke: f!(du, u0, p, t) must execute without error
    du = similar(prob.u0)
    @test_nowarn prob.f(du, prob.u0, prob.p, 0.0)
    @test all(isfinite, du)
end

@testitem "build_ode_problem Path-B: CubedSphere solid-body rotation returns to IC" begin
    using EarthSciDiscretizations: build_ode_problem, CubedSphereGrid, cell_centers
    import SciMLBase
    using OrdinaryDiffEqDefault: solve

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "cubed_sphere_transport.esm")
    gdd_path  = joinpath(repo_root, "test", "fixtures", "curvilinear",
                         "cubed_sphere_transport.gdd.json")

    # Build the spatial operator via Path B.
    prob_template, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob_template isa SciMLBase.ODEProblem
    # GDD: Nc=8 → n_cells = 6 * 8 * 8 = 384
    @test length(prob_template.u0) == 384

    # Inject cosine-bell IC at physical lon/lat coordinates.
    Nc = 8; R = 1.0; h0 = 1.0; r0 = R / 3
    lon0 = 0.0; lat0 = 0.0

    grid = CubedSphereGrid(Nc; R = R)
    lons = cell_centers(grid, :lon)
    lats = cell_centers(grid, :lat)
    N    = 6 * Nc * Nc

    u0_cb = map(1:N) do c
        cosd = sin(lat0) * sin(lats[c]) + cos(lat0) * cos(lats[c]) * cos(lons[c] - lon0)
        r    = R * acos(clamp(cosd, -1.0, 1.0))
        r < r0 ? (h0 / 2) * (1 + cos(π * r / r0)) : 0.0
    end

    # Verify that the RHS is finite and non-zero for the cosine-bell IC.
    du_cb = similar(u0_cb)
    @test_nowarn prob_template.f(du_cb, u0_cb, prob_template.p, 0.0)
    @test all(isfinite, du_cb)
    @test any(!iszero, du_cb)  # transport should produce non-trivial tendency

    # Solve for one full revolution (tspan = (0, 2π) from the ESM domain).
    # remake preserves MTK's initialization data while swapping the IC vector.
    prob = SciMLBase.remake(prob_template; u0 = u0_cb)
    sol = solve(prob; reltol = 1e-6, abstol = 1e-8, save_everystep = false)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # After one revolution the field should return close to the IC.
    u_final   = sol.u[end]
    linf_norm = maximum(abs.(u_final .- u0_cb)) / max(h0, 1e-12)
    @test linf_norm <= 0.5  # ≤ 50% of h0 on Nc=8; tighten in follow-up
end

@testitem "build_ode_problem Path-B: cartesian GDD still routes Path-A" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n16.gdd.json")

    # Should still work via Path A
    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem
    @test length(prob.u0) == 16
    @test haskey(var_map, "u[1]")
end
