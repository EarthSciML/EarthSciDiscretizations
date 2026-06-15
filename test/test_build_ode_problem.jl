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

@testitem "build_ode_problem: expression IC sin(2π·x) loaded from .esm into prob.u0 (esd-wt7)" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "centered_2nd_deriv_uniform", "fixtures", "integration",
                         "diffusion_1d_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_1d_periodic_n32.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N  = 32
    dx = 1.0 / N
    for i in 1:N
        expected = sin(2π * (i - 0.5) * dx)
        @test prob.u0[var_map["u[$i]"]] ≈ expected  rtol = 1e-15
    end
end

@testitem "build_ode_problem: expression IC sin(2π·x)·sin(2π·y) in 2D prob.u0 (esd-wt7)" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "laplacian_2nd_uniform_cartesian", "fixtures", "integration",
                         "diffusion_2d_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "gdd",
                         "cartesian_2d_periodic_n32.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N  = 32
    dx = 1.0 / N
    for i in 1:N, j in 1:N
        expected = sin(2π * (i - 0.5) * dx) * sin(2π * (j - 0.5) * dx)
        @test prob.u0[var_map["u[$i,$j]"]] ≈ expected  rtol = 1e-15
    end
end

@testitem "build_ode_problem: expression IC cos(2π·x) resolution-independent (esd-wt7)" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "upwind_1st", "fixtures", "integration", "advection_1d_pde.esm")

    for (gdd_file, N) in [("adv_n16.gdd.json", 16), ("adv_n32.gdd.json", 32)]
        gdd_path = joinpath(repo_root, "discretizations", "finite_difference",
                            "upwind_1st", "fixtures", "integration", gdd_file)
        prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
        dx = 1.0 / N
        for i in 1:N
            expected = cos(2π * (i - 0.5) * dx)
            @test prob.u0[var_map["u[$i]"]] ≈ expected  rtol = 1e-15
        end
    end
end

# ---------------------------------------------------------------------------
# esd-aij: DUO nearest-neighbor diffusion via pre-bound connectivity
# ---------------------------------------------------------------------------

@testitem "build_ode_problem: DUO nn_diffusion yields finite du on icosahedral level-2 (esd-aij)" begin
    using EarthSciDiscretizations: build_ode_problem
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_duo", "fixtures", "integration",
                         "duo_nn_diffusion_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_duo", "fixtures", "integration",
                         "duo_level2.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    @test prob isa SciMLBase.ODEProblem

    N = 320  # 20 * 4^2 cells at level 2
    @test length(prob.u0) >= N
    @test haskey(var_map, "u[1]")
    @test haskey(var_map, "u[$N]")

    # Evaluate f!(du, u0, p, 0): uniform field → zero Laplacian (du = 0).
    du = similar(prob.u0)
    u0 = copy(prob.u0)
    for c in 1:N
        u0[var_map["u[$c]"]] = 1.0
    end
    prob.f(du, u0, prob.p, 0.0)

    for c in 1:N
        @test isfinite(du[var_map["u[$c]"]])
    end
    for c in 1:N
        @test du[var_map["u[$c]"]] ≈ 0.0  atol = 1e-8
    end
end

@testitem "build_ode_problem: DUO nn_diffusion f! is nonzero for non-uniform IC (esd-aij)" begin
    using EarthSciDiscretizations: build_ode_problem

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_duo", "fixtures", "integration",
                         "duo_nn_diffusion_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_duo", "fixtures", "integration",
                         "duo_level2.gdd.json")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    N  = 320
    du = similar(prob.u0)
    u0 = copy(prob.u0)
    # Non-uniform IC: alternating 0 / 1 by cell parity.
    for c in 1:N
        u0[var_map["u[$c]"]] = iseven(c) ? 1.0 : 0.0
    end
    prob.f(du, u0, prob.p, 0.0)

    # At least some cells should have nonzero tendency.
    @test any(!=(0.0), [du[var_map["u[$c]"]] for c in 1:N])
end

# ---------------------------------------------------------------------------
# esd-aij: MPAS nearest-neighbor diffusion via pre-bound connectivity
# ---------------------------------------------------------------------------

@testitem "build_ode_problem: MPAS nn_diffusion uniform field → zero du (esd-aij)" begin
    using EarthSciDiscretizations: build_ode_problem, mpas_mesh_data
    import SciMLBase

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_mpas", "fixtures", "integration",
                         "mpas_nn_diffusion_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_mpas", "fixtures", "integration",
                         "mpas_triangle_3cell.gdd.json")

    # Synthetic 3-cell triangular mesh.
    # Topology: cells 1,2,3 form a triangle; edges 1=(1,2), 2=(1,3), 3=(2,3).
    # All areas=1, dv_edge=1, dc_edge=1 → Laplacian weight per edge = 1/(1*1).
    # For uniform u=1: du[c] = 2*1 - 2*1 = 0 (each cell has 2 neighbors).
    synth_mesh = mpas_mesh_data(
        lon_cell         = zeros(3),
        lat_cell         = zeros(3),
        area_cell        = [1.0, 1.0, 1.0],
        n_edges_on_cell  = [2, 2, 2],
        cells_on_cell    = [2 1 1; 3 3 2; 0 0 0],  # max_edges × n_cells
        edges_on_cell    = [1 1 2; 2 3 3; 0 0 0],  # max_edges × n_cells
        lon_edge         = zeros(3),
        lat_edge         = zeros(3),
        cells_on_edge    = [1 1 2; 2 3 3],          # 2 × n_edges
        dc_edge          = [1.0, 1.0, 1.0],
        dv_edge          = [1.0, 1.0, 1.0],
        max_edges        = 3,
    )

    prob, var_map = build_ode_problem(esm_path;
                                      grid_ref  = gdd_path,
                                      reader_fn = _ -> synth_mesh)

    @test prob isa SciMLBase.ODEProblem
    @test haskey(var_map, "u[1]") && haskey(var_map, "u[3]")

    du = similar(prob.u0)
    u0 = copy(prob.u0)
    for c in 1:3
        u0[var_map["u[$c]"]] = 1.0
    end
    prob.f(du, u0, prob.p, 0.0)

    # Uniform field must have zero Laplacian.
    for c in 1:3
        @test du[var_map["u[$c]"]] ≈ 0.0  atol = 1e-12
    end
end

@testitem "build_ode_problem: MPAS nn_diffusion Laplacian value matches analytic (esd-aij)" begin
    using EarthSciDiscretizations: build_ode_problem, mpas_mesh_data

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    esm_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_mpas", "fixtures", "integration",
                         "mpas_nn_diffusion_pde.esm")
    gdd_path  = joinpath(repo_root, "discretizations", "finite_difference",
                         "nn_diffusion_mpas", "fixtures", "integration",
                         "mpas_triangle_3cell.gdd.json")

    synth_mesh = mpas_mesh_data(
        lon_cell         = zeros(3),
        lat_cell         = zeros(3),
        area_cell        = [1.0, 1.0, 1.0],
        n_edges_on_cell  = [2, 2, 2],
        cells_on_cell    = [2 1 1; 3 3 2; 0 0 0],
        edges_on_cell    = [1 1 2; 2 3 3; 0 0 0],
        lon_edge         = zeros(3),
        lat_edge         = zeros(3),
        cells_on_edge    = [1 1 2; 2 3 3],
        dc_edge          = [1.0, 1.0, 1.0],
        dv_edge          = [1.0, 1.0, 1.0],
        max_edges        = 3,
    )

    prob, var_map = build_ode_problem(esm_path;
                                      grid_ref  = gdd_path,
                                      reader_fn = _ -> synth_mesh)

    du = similar(prob.u0)
    u0 = copy(prob.u0)
    # IC: u = [0, 1, 0]
    u0[var_map["u[1]"]] = 0.0
    u0[var_map["u[2]"]] = 1.0
    u0[var_map["u[3]"]] = 0.0
    prob.f(du, u0, prob.p, 0.0)

    # Expected Laplacian (dv/dc/area=1 for every edge):
    #   du[1] = u[2] + u[3] - 2*u[1] = 1 + 0 - 0 = 1
    #   du[2] = u[1] + u[3] - 2*u[2] = 0 + 0 - 2 = -2
    #   du[3] = u[1] + u[2] - 2*u[3] = 0 + 1 - 0 = 1
    @test du[var_map["u[1]"]] ≈  1.0  atol = 1e-12
    @test du[var_map["u[2]"]] ≈ -2.0  atol = 1e-12
    @test du[var_map["u[3]"]] ≈  1.0  atol = 1e-12
end
