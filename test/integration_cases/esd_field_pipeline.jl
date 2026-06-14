using EarthSciDiscretizations: build_ode_problem
import SciMLBase
using OrdinaryDiffEqDefault: solve

"""
    run_esd_field_pipeline(name, manifest, base_dir) -> Tuple{Symbol,String}

Dispatch to a sub-runner based on `manifest["pipeline_kind"]`.

Supported pipeline_kinds:
- `diffusion_1d_analytic`   — 1D heat equation, discrete eigenvalue verification
- `diffusion_2d_analytic`   — 2D heat equation, separable discrete eigenvalue
- `advection_1d_convergence` — 1D upwind advection, convergence order check
"""
function run_esd_field_pipeline(name::AbstractString, manifest::AbstractDict,
                                base_dir::AbstractString)
    pkind = String(get(manifest, "pipeline_kind", ""))
    if pkind == "diffusion_1d_analytic"
        return _diffusion_1d_analytic(name, manifest, base_dir)
    elseif pkind == "diffusion_2d_analytic"
        return _diffusion_2d_analytic(name, manifest, base_dir)
    elseif pkind == "advection_1d_convergence"
        return _advection_1d_convergence(name, manifest, base_dir)
    else
        return (:fail, "$name: unknown pipeline_kind '$pkind'")
    end
end

# ---------------------------------------------------------------------------
# 1D diffusion — verify against discrete-eigenvalue exact solution
# ---------------------------------------------------------------------------

function _diffusion_1d_analytic(name::AbstractString, manifest::AbstractDict,
                                base_dir::AbstractString)
    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    gdd_path = abspath(joinpath(base_dir, String(manifest["gdd"])))
    N        = Int(manifest["n_cells"])
    t_final  = Float64(manifest["t_final"])
    tol_max  = Float64(manifest["tolerance"]["max"])

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    dx = 1.0 / N
    u0 = copy(prob.u0)
    for i in 1:N
        u0[var_map["u[$i]"]] = sin(2π * (i - 0.5) * dx)
    end
    prob2 = SciMLBase.ODEProblem(prob.f, u0, (0.0, t_final), prob.p)
    sol = solve(prob2; reltol = 1e-10, abstol = 1e-12, save_everystep = false)

    # Discrete centered-FD second-derivative eigenvalue for sin(2π x):
    #   L_h[sin(2πx)] = 2(cos(2π dx) - 1)/dx² · sin(2πx)
    # The discrete ODE system's exact solution is u_i(t) = sin(2πx_i) exp(λ_h t).
    lambda_h = 2.0 * (cos(2π * dx) - 1.0) / dx^2
    u_final  = sol.u[end]
    linf = maximum(
        abs(u_final[var_map["u[$i]"]] - sin(2π * (i - 0.5) * dx) * exp(lambda_h * t_final))
        for i in 1:N
    )

    if linf <= tol_max
        return (:pass, "$name: L∞ = $(round(linf; sigdigits = 3)) ≤ $tol_max")
    else
        return (:fail, "$name: L∞ = $(round(linf; sigdigits = 3)) exceeds tolerance $tol_max")
    end
end

# ---------------------------------------------------------------------------
# 2D diffusion — separable IC, discrete 5-point Laplacian eigenvalue
# ---------------------------------------------------------------------------

function _diffusion_2d_analytic(name::AbstractString, manifest::AbstractDict,
                                base_dir::AbstractString)
    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    gdd_path = abspath(joinpath(base_dir, String(manifest["gdd"])))
    N        = Int(manifest["n_cells"])
    t_final  = Float64(manifest["t_final"])
    tol_max  = Float64(manifest["tolerance"]["max"])

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    dx = 1.0 / N
    dy = dx  # square grid
    u0 = copy(prob.u0)
    for i in 1:N, j in 1:N
        x_i = (i - 0.5) * dx
        y_j = (j - 0.5) * dy
        u0[var_map["u[$i,$j]"]] = sin(2π * x_i) * sin(2π * y_j)
    end
    prob2 = SciMLBase.ODEProblem(prob.f, u0, (0.0, t_final), prob.p)
    sol = solve(prob2; reltol = 1e-10, abstol = 1e-12, save_everystep = false)

    # Discrete 5-point Laplacian eigenvalue for sin(2πx)sin(2πy) on square grid:
    #   λ_h,2d = 2 · 2(cos(2π dx) - 1)/dx²  (separable, dx = dy)
    lambda_h_1d = 2.0 * (cos(2π * dx) - 1.0) / dx^2
    lambda_h_2d = 2.0 * lambda_h_1d
    u_final = sol.u[end]
    linf = maximum(
        abs(u_final[var_map["u[$i,$j]"]] -
            sin(2π * (i - 0.5) * dx) * sin(2π * (j - 0.5) * dy) * exp(lambda_h_2d * t_final))
        for i in 1:N, j in 1:N
    )

    if linf <= tol_max
        return (:pass, "$name: L∞ = $(round(linf; sigdigits = 3)) ≤ $tol_max")
    else
        return (:fail, "$name: L∞ = $(round(linf; sigdigits = 3)) exceeds tolerance $tol_max")
    end
end

# ---------------------------------------------------------------------------
# 1D advection convergence sweep — first-order upwind, O(h) check
# ---------------------------------------------------------------------------

function _advection_1d_convergence(name::AbstractString, manifest::AbstractDict,
                                   base_dir::AbstractString)
    esm_path  = abspath(joinpath(base_dir, String(manifest["esm"])))
    grid_refs = String.(manifest["grid_refs"])
    v         = Float64(get(manifest, "velocity", 1.0))
    t_final   = Float64(manifest["t_final"])
    min_order = Float64(manifest["convergence"]["min_order"])

    if length(grid_refs) < 2
        return (:fail, "$name: need at least 2 grid_refs for convergence check (got $(length(grid_refs)))")
    end

    errors = Float64[]
    ns     = Int[]
    for gdd_rel in grid_refs
        gdd_path = abspath(joinpath(base_dir, gdd_rel))
        prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

        N  = length(prob.u0)
        dx = 1.0 / N
        u0 = copy(prob.u0)
        for i in 1:N
            u0[var_map["u[$i]"]] = cos(2π * (i - 0.5) * dx)
        end
        prob2 = SciMLBase.ODEProblem(prob.f, u0, (0.0, t_final), prob.p)
        sol = solve(prob2; reltol = 1e-10, abstol = 1e-12, save_everystep = false)

        # Exact: cos(2π(x - v·t)) — periodic translation of the IC
        u_final = sol.u[end]
        linf = maximum(
            abs(u_final[var_map["u[$i]"]] - cos(2π * ((i - 0.5) * dx - v * t_final)))
            for i in 1:N
        )
        push!(errors, linf)
        push!(ns, N)
    end

    # grid_refs ordered coarse→fine; check log2(err_coarse / err_fine) ≥ min_order
    order = log2(errors[1] / errors[2])

    if order >= min_order
        return (:pass,
                "$name: order = $(round(order; sigdigits = 3)) ≥ $min_order " *
                "(N=$(ns[1])→$(ns[2]), L∞=$(round.(errors; sigdigits=3)))")
    else
        return (:fail,
                "$name: order = $(round(order; sigdigits = 3)) < $min_order " *
                "(N=$(ns[1])→$(ns[2]), L∞=$(round.(errors; sigdigits=3)))")
    end
end
