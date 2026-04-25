using EarthSciSerialization: discretize, build_evaluator
using OrdinaryDiffEqDefault: ODEProblem, solve

"""
    run_cartesian_full_pipeline(name, manifest) -> Tuple{Symbol,String}

End-to-end Cartesian benchmark: build a per-cell-scalar ESM model from the
declarative case manifest (`gaussian_advection.esm`), then run the full ESS
pipeline:

    in-memory ESM model
        → EarthSciSerialization.discretize       (rule engine + canonicalization + DAE check)
        → EarthSciSerialization.build_evaluator  (AST → f!(du, u, p, t))
        → OrdinaryDiffEqDefault.solve            (ODE integration)
        → L∞ error vs the analytic translated initial condition

Manifest fields consumed:
- `grid.n`, `grid.extent`, `grid.periodic` (1D periodic uniform)
- `initial_condition.{amplitude, sigma, x0}` (Gaussian)
- `velocity` (scalar, m/s)
- `time.t_final`
- `tolerance.max` (L∞ ceiling)
"""
function run_cartesian_full_pipeline(name::AbstractString, manifest::AbstractDict)
    grid = manifest["grid"]
    if Int(get(grid, "ndim", 0)) != 1
        return (:fail, "$name: only 1D Cartesian supported in this scaffold")
    end
    Nx = Int(grid["n"][1])
    extent = grid["extent"][1]
    L = Float64(extent[2]) - Float64(extent[1])
    dx = L / Nx

    ic = manifest["initial_condition"]
    amp = Float64(ic["amplitude"])
    sigma = Float64(ic["sigma"])
    x0 = Float64(ic["x0"])

    v = Float64(manifest["velocity"])
    t_final = Float64(manifest["time"]["t_final"])
    tol_max = Float64(manifest["tolerance"]["max"])

    cell_x(i) = (i - 0.5) * dx
    gaussian(x) = amp * exp(-((x - x0)^2) / (2 * sigma^2))

    # Build a per-cell-scalar ESM model. Wrap indices so ESS sees periodic BC
    # via the variable name table, not via a `bc` op (which the tree-walk
    # evaluator does not handle). The rule-engine pass through `discretize`
    # canonicalizes equations even when no rules rewrite them.
    wrap(i) = i < 1 ? Nx : (i > Nx ? 1 : i)
    var_name(i) = "u_$(i)"

    variables = Dict{String, Any}()
    for i in 1:Nx
        variables[var_name(i)] = Dict{String, Any}(
            "type" => "state",
            "default" => gaussian(cell_x(i)),
            "units" => "1",
        )
    end
    variables["v"] = Dict{String, Any}(
        "type" => "parameter", "default" => v, "units" => "1",
    )

    equations = Any[]
    for i in 1:Nx
        # Upwind first-order (v > 0): du/dt = -v * (u_i - u_{i-1}) / dx.
        # Use literal -1 multiplier rather than unary `-` so canonicalization
        # does not lower it to the `neg` op (tree-walker has no handler).
        upstream = var_name(wrap(i - 1))
        rhs = Dict{String, Any}(
            "op" => "/", "args" => Any[
                Dict{String, Any}(
                    "op" => "*", "args" => Any[
                        -1,
                        "v",
                        Dict{String, Any}(
                            "op" => "-", "args" => Any[var_name(i), upstream],
                        ),
                    ],
                ),
                dx,
            ],
        )
        push!(
            equations, Dict{String, Any}(
                "lhs" => Dict{String, Any}("op" => "D", "args" => Any[var_name(i)], "wrt" => "t"),
                "rhs" => rhs,
            )
        )
    end

    esm = Dict{String, Any}(
        "esm" => "0.2.0",
        "metadata" => Dict{String, Any}("name" => "gaussian_advection_cartesian_1d"),
        "models" => Dict{String, Any}(
            "M" => Dict{String, Any}(
                "variables" => variables,
                "equations" => equations,
            ),
        ),
    )

    # ESS pipeline: rule engine + AST + solver.
    discretized = discretize(esm)
    f!, u0, p, _tspan, var_map = build_evaluator(discretized)

    prob = ODEProblem(f!, u0, (0.0, t_final), p)
    sol = solve(prob; reltol = 1.0e-6, abstol = 1.0e-8, save_everystep = false)

    u_final = sol.u[end]
    # After one full revolution at v*t_final = L, the analytic solution returns
    # to the IC. (For non-revolution times, shift by v*t modulo L.)
    shift = v * t_final
    analytic = zeros(Nx)
    for i in 1:Nx
        x = cell_x(i)
        x_src = mod(x - shift, L) + Float64(extent[1])
        analytic[i] = gaussian(x_src)
    end
    linf = maximum(abs.(u_final[var_map[var_name(i)]] - analytic[i]) for i in 1:Nx)

    if linf <= tol_max
        return (:pass, "$name: L∞ = $(round(linf; sigdigits = 3)) ≤ $(tol_max) over t=$(t_final)")
    else
        return (:fail, "$name: L∞ = $(round(linf; sigdigits = 3)) exceeds tolerance $(tol_max)")
    end
end
