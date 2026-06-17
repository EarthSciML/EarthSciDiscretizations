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
    esm_path = joinpath(
        repo_root, "test", "fixtures", "curvilinear",
        "latlon_diffusion.esm"
    )
    gdd_path = joinpath(
        repo_root, "test", "fixtures", "curvilinear",
        "latlon_diffusion.gdd.json"
    )

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

@testitem "build_ode_problem Path-B: cartesian GDD still routes Path-A" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd",
        "cartesian_1d_periodic_n16.gdd.json"
    )

    # Should still work via Path A
    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem
    @test length(prob.u0) == 16
    @test haskey(var_map, "u[1]")
end
