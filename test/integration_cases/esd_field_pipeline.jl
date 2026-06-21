using EarthSciDiscretizations: build_ode_problem, build_duo_grid, build_mpas_grid
import SciMLBase
import JSON
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
function run_esd_field_pipeline(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
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
# Returns a vector of NamedTuples: (name, lo, hi, h, N, cell_widths).
# For uniform axes (grid_spacing key): cell_widths is nothing.
# For non-uniform axes (levels key): cell_widths is a Float64 vector of per-cell widths.
# Axes are sorted by name for deterministic ordering.
function _read_gdd_axes(gdd_path::AbstractString)
    gdd = JSON.parse(read(gdd_path, String))
    for (_, domain_spec) in get(gdd, "grids", Dict{String, Any}())
        spatial = get(domain_spec, "spatial", Dict{String, Any}())
        sorted_names = sort!([String(k) for k in keys(spatial)])
        return map(sorted_names) do axis_name
            s = spatial[axis_name]
            if haskey(s, "levels")
                levels = Float64.(s["levels"])
                N = length(levels) - 1
                widths = [levels[k + 1] - levels[k] for k in 1:N]
                lo = levels[1]; hi = levels[end]
                h = (hi - lo) / N
                (
                    name = axis_name, lo = lo, hi = hi, h = h, N = N,
                    cell_widths = widths,
                )
            else
                lo = Float64(get(s, "min", 0.0))
                hi = Float64(get(s, "max", 1.0))
                h = Float64(s["grid_spacing"])
                N = round(Int, (hi - lo) / h)
                (
                    name = axis_name, lo = lo, hi = hi, h = h, N = N,
                    cell_widths = nothing,
                )
            end
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
        v = Float64(get(manifest, "velocity", 1.0))
        maximum(
            let x = ax.lo + (i - 0.5) * ax.h
                    abs(u_vec[var_map["u[$i]"]] - cos(2π * (x - v * t)))
            end
                for i in 1:ax.N
        )
    end,

    # 1-D vertical diffusion with zero Dirichlet BCs at top and bottom.
    # IC: sin(π z).  Exact continuous solution: sin(π z) · exp(−π² t).
    # Works for both uniform and non-uniform grids (reads cell centres from axes).
    # Cell centres: for uniform grids, lo + (k-0.5)*h; for non-uniform,
    # derived from cumulative sum of cell_widths.
    "sin_pi_z_dirichlet" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]
        centres = if ax.cell_widths !== nothing
            ws = ax.cell_widths
            cumws = [sum(ws[1:(k - 1)]) + ws[k] / 2.0 for k in 1:length(ws)]
            cumws
        else
            [ax.lo + (k - 0.5) * ax.h for k in 1:ax.N]
        end
        factor = exp(-π^2 * t)
        maximum(
            abs(u_vec[var_map["u[$k]"]] - sin(π * centres[k]) * factor)
                for k in 1:ax.N
        )
    end,

    # 1-D vertical diffusion, interior cells only.
    # Same PDE and analytic reference as sin_pi_z_dirichlet, but skips
    # `boundary_skip` cells (default 3) at each end.  This isolates the
    # O(h²) interior truncation error from the O(1/h) ghost-cell forcing
    # error that zero-Dirichlet ghost cells introduce at k=1 and k=N.
    # The ghost-cell contamination decays as O(t^k / h^{2k-1}) at cell k
    # from the boundary; with boundary_skip=3 and t_final≤1e-5 it is
    # negligible vs. the O(h²)·t interior truncation error for N≥16.
    "sin_pi_z_dirichlet_interior" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]
        N = ax.N
        n_skip = Int(get(manifest, "boundary_skip", 3))
        k_lo = n_skip + 1
        k_hi = N - n_skip
        k_hi >= k_lo ||
            error("sin_pi_z_dirichlet_interior: no interior cells (N=$N, boundary_skip=$n_skip)")
        centres = if ax.cell_widths !== nothing
            ws = ax.cell_widths
            cumws = cumsum(ws)
            [cumws[k] - ws[k] / 2.0 for k in 1:N]
        else
            [ax.lo + (k - 0.5) * ax.h for k in 1:N]
        end
        factor = exp(-π^2 * t)
        maximum(
            abs(u_vec[var_map["u[$k]"]] - sin(π * centres[k]) * factor)
                for k in k_lo:k_hi
        )
    end,

    # 2-D lat-lon Y_{2,0} colatitude advection (unit sphere, R=1).
    # ESM shape ["lat","lon"] → var_map key "u[i_lat,j_lon]".
    # Axes sorted alphabetically: axes[1]=lat, axes[2]=lon.
    # IC: Y_{2,0}(lat) = 0.25*sqrt(5/pi)*(3*cos^2(lat)-1), lat in [0,pi] colatitude.
    # Exact solution: Y_{2,0}(lat + t).  Pole rows i=1 and i=N_lat skipped
    # (interior-lat scope per bead esd-sqb); lon terms cancel for this lon-independent IC.
    "Y_2_0_latlon_advection_lat" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        lat_ax = axes[findfirst(a -> a.name == "lat", axes)]
        lon_ax = axes[findfirst(a -> a.name == "lon", axes)]
        N_lat = lat_ax.N
        N_lon = lon_ax.N
        c = 0.25 * sqrt(5.0 / pi)
        err = 0.0
        for i in 2:(N_lat - 1)
            lat_c = lat_ax.lo + (i - 0.5) * lat_ax.h
            exact = c * (3 * cos(lat_c + t)^2 - 1)
            for j in 1:N_lon
                idx = get(var_map, "u[$i,$j]", nothing)
                idx === nothing && continue
                err = max(err, abs(u_vec[idx] - exact))
            end
        end
        err
    end,

    # 2-D covariant-FV Laplace-Beltrami diffusion on a regular global lon-lat
    # sphere (Path A, esd-6g4.10). ESM shape ["lat","lon"] → var_map key
    # "u[i_lat,j_lon]"; axes sorted alphabetically: axes[1]=lat, axes[2]=lon.
    # IC: the degree-1 eigenfunction sin(latitude) (the ESM writes sin(lat-π/2) so
    # the no-offset coord sampling materializes sin(true latitude)); under the
    # discrete Laplace-Beltrami it decays as exp(−2t/R²). Exact: sin(lat)·e^{−2t/R²}
    # with lat = lat_ax.lo + (i−0.5)·h (lo = −π/2 ⇒ true latitude, matching the
    # LatLonGrid metric). Pole rows perturb the interior via the zero-ghost
    # boundary, so `boundary_skip` (default 1) interior lat rows are dropped at
    # each pole; the manifest's small t_final keeps the eigen-decay's spatial
    # truncation O(h²)-dominated (cf. sin_pi_z_dirichlet_interior).
    "sin_lat_latlon_laplacian_interior" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        lat_ax = axes[findfirst(a -> a.name == "lat", axes)]
        lon_ax = axes[findfirst(a -> a.name == "lon", axes)]
        N_lat = lat_ax.N
        N_lon = lon_ax.N
        R = Float64(get(manifest, "R", 1.0))
        n_skip = Int(get(manifest, "boundary_skip", 1))
        factor = exp(-2.0 * t / R^2)
        err = 0.0
        for i in (n_skip + 1):(N_lat - n_skip)
            lat_c = lat_ax.lo + (i - 0.5) * lat_ax.h
            exact = sin(lat_c) * factor
            for j in 1:N_lon
                idx = get(var_map, "u[$i,$j]", nothing)
                idx === nothing && continue
                err = max(err, abs(u_vec[idx] - exact))
            end
        end
        err
    end,

    # 2-D Arakawa C-grid divergence via vec_sincos_2d_periodic MMS.
    # F(x,y) = (sin(2πx)cos(2πy), cos(2πx)sin(2πy))  →  ∇·F = 4π cos(2πx)cos(2πy).
    # State variable `div` starts at 0 and accumulates RHS = ∇·F over time t.
    # Exact solution: div(t)[i,j] = t · 4π cos(2πx_c) cos(2πy_c) at cell centre.
    "vec_sincos_2d_periodic" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]; ay = axes[2]
        maximum(
            let x = ax.lo + (i - 0.5) * ax.h,
                    y = ay.lo + (j - 0.5) * ay.h,
                    exact = t * 4π * cos(2π * x) * cos(2π * y)
                    abs(u_vec[var_map["div[$i,$j]"]] - exact)
            end
                for i in 1:ax.N, j in 1:ay.N
        )
    end,

    # 1-D nonlinear (variable-coefficient) diffusion accumulator.
    # Freeze u=sin(2πx) and f=2+sin(2πx); a state `acc` integrates the operator
    # grad(f·grad(u)), so acc(t) = t·grad(f·grad(u)).  Closed form:
    #   grad(f·grad(u)) = (2π)²·(cos(4πx) − 2·sin(2πx)).
    # Reads the accumulator at cell centres (var_map key "acc[i]").  The
    # linear-in-time accumulation makes the ODE-solver error negligible, so the
    # L∞ residual is the rule's O(h²) spatial truncation.
    "sin2pi_nonlinear_diffusion_accum" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]
        maximum(
            let x = ax.lo + (i - 0.5) * ax.h,
                    exact = t * (2π)^2 * (cos(4π * x) - 2 * sin(2π * x))
                    abs(u_vec[var_map["acc[$i]"]] - exact)
            end
                for i in 1:ax.N
        )
    end,

    # 1-D WENO5 flux-divergence accumulator.
    # Freeze U=1 and q=sin(2πx+1); a state `acc` integrates div(U·q), so
    # acc(t) = t·div(U·q) = t·2π·cos(2πx+1).  Reads the accumulator at cell
    # centres (var_map key "acc[i]").  The linear-in-time accumulation isolates
    # the rule's O(h⁵) spatial truncation from the ODE solver.
    "sin2pi_shift1_div_accum" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]
        maximum(
            let x = ax.lo + (i - 0.5) * ax.h,
                    exact = t * 2π * cos(2π * x + 1.0)
                    abs(u_vec[var_map["acc[$i]"]] - exact)
            end
                for i in 1:ax.N
        )
    end,

    # 2-D mixed-derivative accumulator.
    # Freeze u=cos(2π(x+y)); a state `acc` integrates d²u/dxdy, so
    # acc(t) = t·d²u/dxdy = t·(−(2π)²)·cos(2π(x+y)).  Axes sort alphabetically
    # (x, y); reads the accumulator at cell centres (var_map key "acc[i,j]").
    # Freeze-accumulate avoids integrating the bare (indefinite-spectrum) mixed
    # operator; the L∞ residual is the rule's O(h²) spatial truncation.
    "mixed_deriv_2d_accum" =>
        (u_vec, var_map, axes, t, manifest) -> begin
        ax = axes[1]; ay = axes[2]
        maximum(
            let x = ax.lo + (i - 0.5) * ax.h,
                    y = ay.lo + (j - 0.5) * ay.h,
                    exact = t * (-(2π)^2) * cos(2π * (x + y))
                    abs(u_vec[var_map["acc[$i,$j]"]] - exact)
            end
                for i in 1:ax.N, j in 1:ay.N
        )
    end,
)

# ---------------------------------------------------------------------------
# Unstructured-sphere MMS catalog: mms_kind → function
#   (u_vec, var_map, lon_c, lat_c, R_val, t, manifest) -> L∞ error
#
# Unlike _MMS_CONV_CATALOG (which indexes via Cartesian axes), the
# unstructured catalog receives per-cell lon/lat arrays extracted from the
# loaded grid, plus R_val (sphere radius in metres) for the exact decay factor.
# ---------------------------------------------------------------------------
const _UNSTRUCTURED_MMS_CONV_CATALOG = Dict{String, Function}(
    # Spherical harmonic Y₁₁ MMS: u(lon,lat) = sin(lon)*cos(lat).
    # Exact Laplacian: Δu = -2*u/R² (unit-sphere eigenvalue -2, scaled by 1/R²).
    # Exact solution at time t: u(lon,lat,t) = sin(lon)*cos(lat) * exp(-2*t/R²).
    "sin_lon_cos_lat" =>
        (u_vec, var_map, lon_c, lat_c, R_val, t, manifest) -> begin
        Nc = length(lon_c)
        decay = exp(-2.0 * t / R_val^2)
        maximum(
            let idx = get(var_map, "u[$c]", nothing)
                    idx === nothing ? 0.0 :
                    abs(u_vec[idx] - sin(lon_c[c]) * cos(lat_c[c]) * decay)
            end
                for c in 1:Nc
        )
    end,
)

function _run_unstructured_mms_convergence(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    grid_refs = String.(manifest["grid_refs"])
    t_final = Float64(manifest["t_final"])
    min_order = Float64(manifest["convergence"]["min_order"])
    mms_kind = String(get(manifest, "mms_kind", ""))

    length(grid_refs) < 2 &&
        return (:fail, "$name: need ≥ 2 grid_refs for convergence check (got $(length(grid_refs)))")

    linf_fn = get(_UNSTRUCTURED_MMS_CONV_CATALOG, mms_kind, nothing)
    linf_fn === nothing &&
        return (:fail, "$name: unknown unstructured mms_kind '$mms_kind'")

    errors = Float64[]
    for gdd_rel in grid_refs
        gdd_path = isabspath(gdd_rel) ? gdd_rel : abspath(joinpath(base_dir, gdd_rel))

        gdd = JSON.parse(read(gdd_path, String))
        domain_spec = get(get(gdd, "grids", Dict{String, Any}()), "domain", Dict{String, Any}())
        family = String(get(domain_spec, "family", ""))

        local lon_c::Vector{Float64}
        local lat_c::Vector{Float64}
        local R_val::Float64
        local Nc::Int

        if family == "mpas"
            loader_spec = get(domain_spec, "loader", nothing)
            loader_path = String(get(loader_spec, "path", ""))
            R_sphere = Float64(get(loader_spec, "sphere_radius", 6.371e6))
            grid = build_mpas_grid(
                loader = Dict(
                    "path" => loader_path,
                    "reader" => String(get(loader_spec, "reader", "auto")),
                    "check" => String(get(loader_spec, "check", "strict"))
                ),
                R = R_sphere,
            )
            lon_c = grid.mesh.lon_cell
            lat_c = grid.mesh.lat_cell
            R_val = Float64(grid.R)
            Nc = length(lon_c)
        elseif family == "duo"
            loader_spec = get(domain_spec, "loader", nothing)
            loader_path = String(get(loader_spec, "path", ""))
            R_sphere = Float64(get(loader_spec, "sphere_radius", 6.371e6))
            grid = build_duo_grid(
                loader = (
                    path = loader_path,
                    reader = String(get(loader_spec, "reader", "auto")),
                ),
                R = R_sphere,
            )
            lon_c = grid.cc_lon
            lat_c = grid.cc_lat
            R_val = Float64(grid.R)
            Nc = length(lon_c)
        else
            return (:fail, "$name: unstructured runner does not support grid family '$family'")
        end

        extra_ics = Dict{String, Float64}("u[$c]" => sin(lon_c[c]) * cos(lat_c[c]) for c in 1:Nc)
        prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path, extra_ics = extra_ics)
        prob2 = SciMLBase.ODEProblem(prob.f, prob.u0, (0.0, t_final), prob.p)
        sol = solve(prob2; reltol = 1.0e-10, abstol = 1.0e-12, save_everystep = false)

        push!(errors, linf_fn(sol.u[end], var_map, lon_c, lat_c, R_val, t_final, manifest))
    end

    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    order = minimum(orders)

    if order >= min_order
        return (
            :pass, "$name: min order = $(round(order; sigdigits = 3)) ≥ $min_order " *
                "(orders=$(round.(orders; sigdigits = 3)), L∞=$(round.(errors; sigdigits = 3)))",
        )
    else
        return (
            :fail, "$name: min order = $(round(order; sigdigits = 3)) < $min_order " *
                "(orders=$(round.(orders; sigdigits = 3)), L∞=$(round.(errors; sigdigits = 3)))",
        )
    end
end

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
function run_eigenvalue_case(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    gdd_rel = String(manifest["gdd"])
    gdd_path = isabspath(gdd_rel) ? gdd_rel : abspath(joinpath(base_dir, gdd_rel))
    t_final = Float64(manifest["t_final"])
    tol_max = Float64(manifest["tolerance"]["max"])
    ev_kind = String(get(manifest, "eigenvalue_kind", "cartesian_1d_laplacian"))

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
    sol = solve(prob2; reltol = 1.0e-10, abstol = 1.0e-12, save_everystep = false)
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
function run_mms_convergence_case(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    grid_refs = String.(manifest["grid_refs"])

    # Detect unstructured family by probing the first GDD. Unstructured grids
    # have no spatial axes (n_cells-only), so dispatch to the dedicated runner.
    if !isempty(grid_refs)
        probe_path = let p = grid_refs[1]
            isabspath(p) ? p : abspath(joinpath(base_dir, p))
        end
        if isfile(probe_path)
            probe_gdd = JSON.parse(read(probe_path, String))
            probe_domain = get(get(probe_gdd, "grids", Dict()), "domain", Dict())
            probe_family = String(get(probe_domain, "family", ""))
            if probe_family in ("mpas", "duo")
                return _run_unstructured_mms_convergence(name, manifest, base_dir)
            end
        end
    end

    esm_path = abspath(joinpath(base_dir, String(manifest["esm"])))
    t_final = Float64(manifest["t_final"])
    min_order = Float64(manifest["convergence"]["min_order"])
    mms_kind = String(get(manifest, "mms_kind", ""))

    if length(grid_refs) < 2
        return (
            :fail,
            "$name: need at least 2 grid_refs for convergence check " *
                "(got $(length(grid_refs)))",
        )
    end

    linf_fn = get(_MMS_CONV_CATALOG, mms_kind, nothing)
    if linf_fn === nothing
        return (:fail, "$name: unknown mms_kind '$mms_kind'")
    end

    errors = Float64[]
    ns = Int[]
    for gdd_rel in grid_refs
        gdd_path = isabspath(gdd_rel) ? gdd_rel : abspath(joinpath(base_dir, gdd_rel))

        axes = _read_gdd_axes(gdd_path)
        if isempty(axes)
            return (:fail, "$name: no spatial axes found in GDD '$gdd_rel'")
        end

        prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)
        prob2 = SciMLBase.ODEProblem(prob.f, prob.u0, (0.0, t_final), prob.p)
        sol = solve(prob2; reltol = 1.0e-10, abstol = 1.0e-12, save_everystep = false)

        push!(errors, linf_fn(sol.u[end], var_map, axes, t_final, manifest))
        push!(ns, axes[1].N)
    end

    order = log2(errors[1] / errors[2])

    if order >= min_order
        return (
            :pass,
            "$name: order = $(round(order; sigdigits = 3)) ≥ $min_order " *
                "(N=$(ns[1])→$(ns[2]), L∞=$(round.(errors; sigdigits = 3)))",
        )
    else
        return (
            :fail,
            "$name: order = $(round(order; sigdigits = 3)) < $min_order " *
                "(N=$(ns[1])→$(ns[2]), L∞=$(round.(errors; sigdigits = 3)))",
        )
    end
end

# ---------------------------------------------------------------------------
# Legacy runners — thin wrappers that delegate to the grid-agnostic helpers.
# Existing manifest schemas are preserved; no callers need updating.
# ---------------------------------------------------------------------------

function _diffusion_1d_analytic(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    return run_eigenvalue_case(
        name,
        merge(Dict{String, Any}("eigenvalue_kind" => "cartesian_1d_laplacian"), manifest),
        base_dir
    )
end

function _diffusion_2d_analytic(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    return run_eigenvalue_case(
        name,
        merge(Dict{String, Any}("eigenvalue_kind" => "cartesian_2d_laplacian"), manifest),
        base_dir
    )
end

function _advection_1d_convergence(
        name::AbstractString, manifest::AbstractDict,
        base_dir::AbstractString
    )
    return run_mms_convergence_case(
        name,
        merge(Dict{String, Any}("mms_kind" => "cos_2pi_x_transport"), manifest),
        base_dir
    )
end
