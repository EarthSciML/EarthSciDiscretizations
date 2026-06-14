using Test
using TestItems
import SciMLBase

@testitem "build_ode_problem: returns ODEProblem and var_map for 1D diffusion" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n16.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem
    @test var_map isa Dict{String, Int}

    # N=16 cells → 16 state entries
    N = 16
    @test length(prob.u0) == N
    @test haskey(var_map, "u[1]")
    @test haskey(var_map, "u[$N]")
end

@testitem "build_ode_problem: f!(du,u0,p,t) reproduces centered-FD Laplacian" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n16.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N  = 16
    dx = 1.0 / N

    # Set a nontrivial initial condition so the check is not trivial.
    u0 = copy(prob.u0)
    for i in 1:N
        u0[var_map["u[$i]"]] = sin(2π * (i - 0.5) * dx)
    end

    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    # The centered_2nd_deriv_uniform stencil produces
    #   du[i] = (u[i-1] - 2*u[i] + u[i+1]) / dx²
    # with periodic wrap.
    for i in 1:N
        im1 = mod1(i - 1, N)
        ip1 = mod1(i + 1, N)
        expected = (u0[var_map["u[$im1]"]] -
                    2 * u0[var_map["u[$i]"]] +
                    u0[var_map["u[$ip1]"]]) / dx^2
        @test du[var_map["u[$i]"]] ≈ expected rtol = 1e-10
    end
end

@testitem "build_ode_problem: no solver invoked, tspan populated" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n16.gdd.json")

    prob, _ = build_ode_problem(esm_path; grid_ref = gdd_path)

    # tspan is a 2-tuple of reals; no solve needed to access it.
    @test prob.tspan isa Tuple{<:Real, <:Real}
    @test prob.tspan[1] < prob.tspan[2]
    # u0 is populated (non-empty)
    @test !isempty(prob.u0)
end

@testitem "build_ode_problem: N=32 grid via different gdd ref" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "test", "fixtures", "diffusion_1d.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n32.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N  = 32
    dx = 1.0 / N

    @test length(prob.u0) == N
    @test haskey(var_map, "u[$N]")

    # Stencil check at N=32
    u0 = copy(prob.u0)
    for i in 1:N
        u0[var_map["u[$i]"]] = Float64(i)
    end
    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    for i in 1:N
        im1 = mod1(i - 1, N)
        ip1 = mod1(i + 1, N)
        expected = (u0[var_map["u[$im1]"]] -
                    2 * u0[var_map["u[$i]"]] +
                    u0[var_map["u[$ip1]"]]) / dx^2
        @test du[var_map["u[$i]"]] ≈ expected rtol = 1e-10
    end
end
