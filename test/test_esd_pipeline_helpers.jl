using Test
using TestItems

@testitem "esd_pipeline: _read_gdd_axes returns correct axes for 1D GDD" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd",
        "cartesian_1d_periodic_n32.gdd.json"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    _read_gdd_axes = Base.invokelatest(getfield, getfield(m, :IntegrationCases), :_read_gdd_axes)

    axes = Base.invokelatest(_read_gdd_axes, gdd_path)

    @test length(axes) == 1
    @test axes[1].name == "x"
    @test axes[1].lo ≈ 0.0
    @test axes[1].hi ≈ 1.0
    @test axes[1].h ≈ 0.03125
    @test axes[1].N == 32
end

@testitem "esd_pipeline: _read_gdd_axes returns correct axes for 2D GDD" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    gdd_path = joinpath(
        repo_root, "discretizations", "gdd",
        "cartesian_2d_periodic_n32.gdd.json"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    _read_gdd_axes = Base.invokelatest(getfield, getfield(m, :IntegrationCases), :_read_gdd_axes)

    axes = Base.invokelatest(_read_gdd_axes, gdd_path)

    @test length(axes) == 2
    @test axes[1].name == "x"
    @test axes[2].name == "y"
    @test axes[1].N == 32
    @test axes[2].N == 32
    @test axes[1].h ≈ axes[2].h
end

@testitem "esd_pipeline: run_eigenvalue_case passes 1D diffusion (eigenvalue_analytic kind)" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "centered_2nd_deriv_uniform",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_eigenvalue_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_eigenvalue_case
    )

    manifest = Dict{String, Any}(
        "pipeline_kind" => "eigenvalue_analytic",
        "esm" => "diffusion_1d_pde.esm",
        "gdd" => "../../../../gdd/cartesian_1d_periodic_n32.gdd.json",
        "eigenvalue_kind" => "cartesian_1d_laplacian",
        "t_final" => 0.01,
        "tolerance" => Dict("metric" => "Linf", "max" => 1.0e-6),
    )

    outcome, msg = Base.invokelatest(run_eigenvalue_case, "test_ev_1d", manifest, base_dir)

    @test outcome === :pass
    @test occursin("L∞", msg)
    @test occursin("≤", msg)
end

@testitem "esd_pipeline: run_eigenvalue_case passes 2D diffusion (eigenvalue_analytic kind)" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "laplacian_2nd_uniform_cartesian",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_eigenvalue_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_eigenvalue_case
    )

    manifest = Dict{String, Any}(
        "pipeline_kind" => "eigenvalue_analytic",
        "esm" => "diffusion_2d_pde.esm",
        "gdd" => "../../../../gdd/cartesian_2d_periodic_n32.gdd.json",
        "eigenvalue_kind" => "cartesian_2d_laplacian",
        "t_final" => 0.005,
        "tolerance" => Dict("metric" => "Linf", "max" => 1.0e-6),
    )

    outcome, msg = Base.invokelatest(run_eigenvalue_case, "test_ev_2d", manifest, base_dir)

    @test outcome === :pass
    @test occursin("L∞", msg)
end

@testitem "esd_pipeline: run_mms_convergence_case passes 1D advection (mms_convergence kind)" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "upwind_1st",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_mms_convergence_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_mms_convergence_case
    )

    manifest = Dict{String, Any}(
        "pipeline_kind" => "mms_convergence",
        "esm" => "advection_1d_pde.esm",
        "grid_refs" => ["adv_n16.gdd.json", "adv_n32.gdd.json"],
        "mms_kind" => "cos_2pi_x_transport",
        "velocity" => 1.0,
        "t_final" => 0.1,
        "convergence" => Dict("min_order" => 0.8),
    )

    outcome, msg = Base.invokelatest(
        run_mms_convergence_case, "test_mms_conv_1d", manifest, base_dir
    )

    @test outcome === :pass
    @test occursin("order", msg)
    @test occursin("≥", msg)
end

@testitem "esd_pipeline: run_eigenvalue_case fails with unknown eigenvalue_kind" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "centered_2nd_deriv_uniform",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_eigenvalue_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_eigenvalue_case
    )

    manifest = Dict{String, Any}(
        "esm" => "diffusion_1d_pde.esm",
        "gdd" => "../../../../gdd/cartesian_1d_periodic_n32.gdd.json",
        "eigenvalue_kind" => "no_such_kind",
        "t_final" => 0.01,
        "tolerance" => Dict("max" => 1.0e-6),
    )

    outcome, msg = Base.invokelatest(run_eigenvalue_case, "bad_ev_kind", manifest, base_dir)

    @test outcome === :fail
    @test occursin("no_such_kind", msg)
end

@testitem "esd_pipeline: run_mms_convergence_case fails with unknown mms_kind" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "upwind_1st",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_mms_convergence_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_mms_convergence_case
    )

    manifest = Dict{String, Any}(
        "esm" => "advection_1d_pde.esm",
        "grid_refs" => ["adv_n16.gdd.json", "adv_n32.gdd.json"],
        "mms_kind" => "no_such_mms",
        "t_final" => 0.1,
        "convergence" => Dict("min_order" => 0.8),
    )

    outcome, msg = Base.invokelatest(
        run_mms_convergence_case, "bad_mms_kind", manifest, base_dir
    )

    @test outcome === :fail
    @test occursin("no_such_mms", msg)
end

@testitem "esd_pipeline: run_mms_convergence_case fails with < 2 grid_refs" begin
    import EarthSciDiscretizations

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    base_dir = joinpath(
        repo_root, "discretizations",
        "finite_difference", "upwind_1st",
        "fixtures", "integration"
    )

    m = Module()
    Base.include(m, joinpath(repo_root, "test", "integration_cases", "IntegrationCases.jl"))
    run_mms_convergence_case = Base.invokelatest(
        getfield, getfield(m, :IntegrationCases), :run_mms_convergence_case
    )

    manifest = Dict{String, Any}(
        "esm" => "advection_1d_pde.esm",
        "grid_refs" => ["adv_n16.gdd.json"],
        "mms_kind" => "cos_2pi_x_transport",
        "t_final" => 0.1,
        "convergence" => Dict("min_order" => 0.8),
    )

    outcome, msg = Base.invokelatest(
        run_mms_convergence_case, "one_grid_ref", manifest, base_dir
    )

    @test outcome === :fail
    @test occursin("2 grid_refs", msg)
end
