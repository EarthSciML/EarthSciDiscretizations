using EarthSciDiscretizations: build_ode_problem
import SciMLBase
using OrdinaryDiffEqDefault: solve

"""
    run_esd_field_pipeline(name, manifest, base_dir) -> Tuple{Symbol,String}

Dispatch to a sub-runner based on `manifest["pipeline_kind"]`.

Supported pipeline_kinds:
- `diffusion_1d_analytic`   — 1D heat equation, discrete eigenvalue verification
  (legacy wrapper; delegates to `run_eigenvalue_case`)
- `diffusion_2d_analytic`   — 2D heat equation, separable discrete eigenvalue
  (legacy wrapper; delegates to `run_eigenvalue_case`)
- `advection_1d_convergence` — 1D upwind advection, convergence order check
  (legacy wrapper; delegates to `run_mms_convergence_case`)
- `eigenvalue_analytic`     — generic discrete-eigenvalue verifier
- `mms_convergence`         — generic MMS convergence sweep
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
    elseif pkind == "eigenvalue_analytic"
        return run_eigenvalue_case(name, manifest, base_dir)
    elseif pkind == "mms_convergence"
        return run_mms_convergence_case(name, manifest, base_dir)
    else
        return (:fail, "$name: unknown pipeline_kind '$pkind'")
    end
end

# ---------------------------------------------------------------------------
# GDD axis helper
# ---------------------------------------------------------------------------

# Read spatial axes from a GDD file (first domain found).
# Returns a vector of NamedTuples: (name, lo, hi, h, N).
# Axes are sorted by name for deterministic ordering.
function _read_gdd_axes(gdd_path::AbstractString)
    gdd = JSON.parse(read(gdd_path, String))
    for (_, domain_spec) in get(gdd, "grids", Dict{String,Any}())
        spatial = get(domain_spec, "spatial", Dict{String,Any}())
        sorted_names = sort!([String(k) for k in keys(spatial)])
        return map(sorted_names) do axis_name
            s  = spatial[axis_name]
            lo = Float64(get(s, "min", 0.0))
            hi = Float64(get(s, "max", 1.0))
            h  = Float64(s["grid_spacing"])
            N  = round(Int, (hi - lo) / h)
            (name = axis_name, lo = lo, hi = hi, h = h, N = N)
        end
    end
    return NamedTuple[]
end

# ---------------------------------------------------------------------------
# Discrete eigenvalue catalog, keyed by eigenvalue_kind.
# Each entry is a function of the GDD axes tuple → Float64 eigenvalue.
# The paired ESM expression IC must be an eigenfunction of the discrete operator.
# ---------------------------------------------------------------------------
const _EV_KIND_CATALOG = Dict{String, Function}(
    # 2nd-order centered FD d²/dx² on a 1-D uniform periodic grid.
    # Eigenvector: sin(2π x).  λ_h = 2(cos(2π h) − 1) / h²
    "cartesian_1d_laplacian" =>
        axes -> (h = axes[1].h; 2.0 * (cos(2π * h) - 1.0) / h^2),

    # 5-point Cartesian Laplacian on a 2-D uniform periodic grid.
    # Eigenvector: sin(2πx)·sin(2πy).  λ_h = λ_x + λ_y (Kronecker sum, separable)
    "cartesian_2d_laplacian" =>
        axes -> begin
            hx = axes[1].h; hy = axes[2].h
            (2.0 * (cos(2π * hx) - 1.0) / hx^2) + (2.0 * (cos(2π * hy) - 1.0) / hy^2)
        end,
)

# ---------------------------------------------------------------------------
# MMS convergence catalog, keyed by mms_kind.
# Each entry: (u_vec, var_map, axes, t_final, manifest) → Float64 L∞ error.
# Cell-centre coordinates are derived from axes (not hardcoded to Cartesian),
# so new grid families extend this catalog without touching the runner.
# ---------------------------------------------------------------------------
const _MMS_CONV_CATALOG = Dict{String, Function}(
    # 1-D periodic transport at constant velocity v.
    # IC: cos(2π x).  Exact: cos(2π(x − v·t)).
    # Manifest field: velocity (Float64, default 1.0).
    "cos_2pi_x_transport" =>
        (u_vec, var_map, axes, t, manifest) -> begin
            ax = axes[1]
            v  = Float64(get(manifest, "velocity", 1.0))
            maximum(
                let x = ax.lo + (i - 0.5) * ax.h
                    abs(u_vec[var_map["u[$i]"]] - cos(2π * (x - v * t)))
                end
                for i in 1:ax.N
            )
        end,
)

# ---------------------------------------------------------------------------
# run_eigenvalue_case — grid-agnostic discrete-eigenvalue verifier
#
# Verifies that the ODE solver error is within tolerance when the IC is an
# eigenfunction of the discrete spatial operator, so the exact ODE solution
# is u(t) = u(0) · exp(λ_h · t).
#
# Grid coordinates are read from the GDD; the IC comes from the ESM's
# expression IC block (evaluated by build_ode_problem). No Cartesian
# assumptions are made in this function.
#
# Manifest fields:
#   esm              path to PDE ESM (relative to base_dir)
#   gdd              path to GDD (relative to base_dir, or absolute)
#   eigenvalue_kind  key into _EV_KIND_CATALOG (default "cartesian_1d_laplacian")
#   t_final          final integration time
#   tolerance.max    maximum acceptable L∞ error
# ---------------------------------------------------------------------------
function run_eigenvalue_case(name::AbstractString, manifest::AbstractDict,
                              base_dir::AbstractString)
    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    gdd_rel  = String(manifest["gdd"])
    gdd_path = isabspath(gdd_rel) ? gdd_rel : abspath(joinpath(base_dir, gdd_rel))
    t_final  = Float64(manifest["t_final"])
    tol_max  = Float64(manifest["tolerance"]["max"])
    ev_kind  = String(get(manifest, "eigenvalue_kind", "cartesian_1d_laplacian"))

    ev_fn = get(_EV_KIND_CATALOG, ev_kind, nothing)
    if ev_fn === nothing
        return (:fail, "$name: unknown eigenvalue_kind '$ev_kind'")
    end

    axes = _read_gdd_axes(gdd_path)
    if isempty(axes)
        return (:fail, "$name: no spatial axes found in GDD '$gdd_rel'")
    end
    lambda_h = ev_fn(axes)

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    prob2 = SciMLBase.ODEProblem(prob.f, prob.u0, (0.0, t_final), prob.p)
    sol = solve(prob2; reltol = 1e-10, abstol = 1e-12, save_everystep = false)
    u_final = sol.u[end]

    # Exact ODE solution: u_i(t) = u_i(0) · exp(λ_h · t) for every cell i.
    factor = exp(lambda_h * t_final)
    linf = maximum(abs(u_final[idx] - prob.u0[idx] * factor) for (_, idx) in var_map)

    if linf <= tol_max
        return (:pass, "$name: L∞ = $(round(linf; sigdigits = 3)) ≤ $tol_max")
    else
        return (:fail, "$name: L∞ = $(round(linf; sigdigits = 3)) exceeds tolerance $tol_max")
    end
end

# ---------------------------------------------------------------------------
# run_mms_convergence_case — grid-agnostic MMS convergence sweep
#
# Runs build_ode_problem for each GDD in grid_refs (coarse → fine), computes
# L∞ error against the analytic reference from the MMS catalog, and checks
# that log₂(err_coarse / err_fine) ≥ min_order.
#
# Cell-centre coordinates are read from each GDD so the helper works for any
# grid family whose GDD declares spatial axes with grid_spacing, min, max.
#
# Manifest fields:
#   esm              path to PDE ESM (relative to base_dir)
#   grid_refs        list of GDD paths, coarse-to-fine (relative to base_dir)
#   mms_kind         key into _MMS_CONV_CATALOG
#   t_final          final integration time
#   convergence.min_order  minimum acceptable convergence order (log₂ ratio)
#   (additional fields may be required by the chosen mms_kind, e.g. velocity)
# ---------------------------------------------------------------------------
function run_mms_convergence_case(name::AbstractString, manifest::AbstractDict,
                                   base_dir::AbstractString)
    esm_path  = abspath(joinpath(base_dir, String(manifest["esm"])))
    grid_refs = String.(manifest["grid_refs"])
    t_final   = Float64(manifest["t_final"])
    min_order = Float64(manifest["convergence"]["min_order"])
    mms_kind  = String(get(manifest, "mms_kind", ""))

    if length(grid_refs) < 2
        return (:fail,
                "$name: need at least 2 grid_refs for convergence check " *
                "(got $(length(grid_refs)))")
    end

    linf_fn = get(_MMS_CONV_CATALOG, mms_kind, nothing)
    if linf_fn === nothing
        return (:fail, "$name: unknown mms_kind '$mms_kind'")
    end

    errors = Float64[]
    ns     = Int[]
    for gdd_rel in grid_refs
        gdd_path = isabspath(gdd_rel) ? gdd_rel : abspath(joinpath(base_dir, gdd_rel))

        axes = _read_gdd_axes(gdd_path)
        if isempty(axes)
            return (:fail, "$name: no spatial axes found in GDD '$gdd_rel'")
        end

        prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
        prob2 = SciMLBase.ODEProblem(prob.f, prob.u0, (0.0, t_final), prob.p)
        sol = solve(prob2; reltol = 1e-10, abstol = 1e-12, save_everystep = false)

        push!(errors, linf_fn(sol.u[end], var_map, axes, t_final, manifest))
        push!(ns, axes[1].N)
    end

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

# ---------------------------------------------------------------------------
# Legacy runners — thin wrappers that delegate to the grid-agnostic helpers.
# Existing manifest schemas are preserved; no callers need updating.
# ---------------------------------------------------------------------------

function _diffusion_1d_analytic(name::AbstractString, manifest::AbstractDict,
                                base_dir::AbstractString)
    return run_eigenvalue_case(name,
        merge(Dict{String,Any}("eigenvalue_kind" => "cartesian_1d_laplacian"), manifest),
        base_dir)
end

function _diffusion_2d_analytic(name::AbstractString, manifest::AbstractDict,
                                base_dir::AbstractString)
    return run_eigenvalue_case(name,
        merge(Dict{String,Any}("eigenvalue_kind" => "cartesian_2d_laplacian"), manifest),
        base_dir)
end

function _advection_1d_convergence(name::AbstractString, manifest::AbstractDict,
                                   base_dir::AbstractString)
    return run_mms_convergence_case(name,
        merge(Dict{String,Any}("mms_kind" => "cos_2pi_x_transport"), manifest),
        base_dir)
end
