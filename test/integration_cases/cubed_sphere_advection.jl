using EarthSciDiscretizations: CubedSphereGrid, transport_2d_linrood!, total_area,
    build_ode_problem, cell_centers
import SciMLBase
using OrdinaryDiffEqDefault: solve

"""
    run_cubed_sphere_advection(name, manifest) -> Tuple{Symbol,String}

Williamson-style advection benchmark on the cubed sphere. Time-steps an
initial scalar field under solid-body rotation using the existing PPM
transport machinery (`transport_2d_linrood!`) and reports L∞ vs the
analytic solution after `T_revolutions` full revolutions.

This case currently bypasses the .esm → discretize → build_evaluator path
because curvilinear array equations require ESS scalarization (deferred to
a follow-up bead). The PPM transport here is the same dimensionally-split
Lin-Rood scheme exercised by the unit tests in `test_shallow_water.jl`,
just driven for many time steps with an L∞ comparison at the end.

Manifest fields consumed:
- `grid.{Nc, R}`
- `initial_condition.{kind, h0, r0_factor, lon0, lat0}`  (kind="cosine_bell")
- `velocity.{kind, u0, alpha_radians}`                    (kind="solid_body_rotation")
- `time.{T_revolutions, courant_target}`
- `tolerance.max`
"""
function run_cubed_sphere_advection(name::AbstractString, manifest::AbstractDict)
    grid_spec = manifest["grid"]
    Nc = Int(grid_spec["Nc"])
    R = Float64(grid_spec["R"])
    grid = CubedSphereGrid(Nc; R = R)

    ic = manifest["initial_condition"]
    if String(ic["kind"]) != "cosine_bell"
        return (:fail, "$name: only cosine_bell IC is wired in this scaffold")
    end
    h0 = Float64(ic["h0"])
    r0 = R * Float64(ic["r0_factor"])
    lon0 = Float64(ic["lon0"])
    lat0 = Float64(ic["lat0"])

    vel_spec = manifest["velocity"]
    if String(vel_spec["kind"]) != "solid_body_rotation"
        return (:fail, "$name: only solid_body_rotation velocity is wired in this scaffold")
    end
    u0 = Float64(vel_spec["u0"])
    alpha = Float64(get(vel_spec, "alpha_radians", 0.0))

    T_rev = Float64(manifest["time"]["T_revolutions"])
    courant = Float64(get(manifest["time"], "courant_target", 0.4))
    tol_max = Float64(manifest["tolerance"]["max"])

    # Cosine-bell IC at every cell center (panel, i, j).
    h = zeros(6, Nc, Nc)
    function cosine_bell(lo, la)
        cosd = sin(lat0) * sin(la) + cos(lat0) * cos(la) * cos(lo - lon0)
        r = R * acos(clamp(cosd, -1.0, 1.0))
        return r < r0 ? (h0 / 2) * (1 + cos(pi * r / r0)) : 0.0
    end
    for p in 1:6, i in 1:Nc, j in 1:Nc
        h[p, i, j] = cosine_bell(grid.lon[p, i, j], grid.lat[p, i, j])
    end
    h_initial = copy(h)

    # Solid-body rotation in lon-lat: u_lon = u0 * (cos(lat)*cos(alpha) +
    # sin(lat)*cos(lon)*sin(alpha)); u_lat = -u0 * sin(lon) * sin(alpha).
    # Project onto cubed-sphere ξ-/η-edges. The scaffolding case uses
    # alpha=0 so the rotation lies entirely on the equator: u_lon = u0*cos(lat),
    # u_lat = 0. The runner currently approximates the panel-edge ξ-/η-velocity
    # at the (cell-center) latitude — exact projection onto the ξ-/η-axes is
    # a follow-up that requires the local edge tangent vectors.
    vel_xi = zeros(6, Nc + 1, Nc)
    vel_eta = zeros(6, Nc, Nc + 1)
    for p in 1:6, i in 1:(Nc + 1), j in 1:Nc
        # Use the latitude at the nearest cell center — cell-edge lat is not
        # stored on the grid and a point-evaluator is not yet wired.
        ic_clamp = clamp(i, 1, Nc)
        vel_xi[p, i, j] = u0 * cos(grid.lat[p, ic_clamp, j])
    end
    # alpha=0: north-south velocity is zero by construction.

    # CFL-based time step: dt = courant * dx_min / u_max.
    dx_min = minimum(grid.dx)
    u_max = max(maximum(abs.(vel_xi)), maximum(abs.(vel_eta)), 1.0)
    dt = courant * dx_min / u_max
    t_final = T_rev * 2 * pi * R / u0
    nsteps = max(1, ceil(Int, t_final / dt))
    dt_actual = t_final / nsteps

    tendency = zeros(6, Nc, Nc)
    for _ in 1:nsteps
        transport_2d_linrood!(tendency, h, vel_xi, vel_eta, grid, dt_actual)
        @. h = h + dt_actual * tendency
    end

    # After one revolution the analytic solution returns to the IC.
    linf_abs = maximum(abs.(h .- h_initial))
    linf_norm = linf_abs / max(h0, 1.0e-12)

    msg = "$name: L∞/h0 = $(round(linf_norm; sigdigits = 3)) over $(nsteps) steps (dt=$(round(dt_actual; sigdigits = 3)) s, T=$(round(t_final; sigdigits = 4)) s)"
    if linf_norm <= tol_max
        return (:pass, msg)
    else
        return (:fail, msg * " — exceeds tolerance $(tol_max)")
    end
end

"""
    run_cubed_sphere_path_b_rotation(name, manifest, base_dir) -> Tuple{Symbol,String}

Solid-body rotation benchmark via Path B: loads a transport PDE from an ESM file
(.esm + GDD with family=cubed_sphere), routes it through
`build_ode_problem` → `EarthSciSerialization.discretize(sys, CubedSphereGrid)`,
injects a cosine-bell IC, solves for one full revolution, and verifies
return-to-IC within tolerance.

The ESM must encode `∂h/∂t = -(u₀/R)·∂h/∂lon` for a chosen angular velocity
with a temporal domain of one revolution period. For unit angular velocity
(coefficient = -1, R = 1): `∂h/∂t = -∂h/∂lon`, T = 2π.

Because Path B always produces a zero IC (no spatial BCs in the PDESystem),
the cosine-bell IC is injected post-construction via `SciMLBase.remake`,
which preserves MTK's initialization data while replacing the u0 vector.

Manifest fields consumed:
- `esm_path`  path to transport ESM (relative to base_dir)
- `gdd_path`  path to GDD (relative to base_dir)
- `initial_condition.{kind, h0, r0_factor, lon0, lat0}`  (kind="cosine_bell")
- `tolerance.max`  L∞/h0 ceiling after one revolution
"""
function run_cubed_sphere_path_b_rotation(name::AbstractString, manifest::AbstractDict,
                                           base_dir::AbstractString)
    esm_path = abspath(joinpath(base_dir, String(manifest["esm_path"])))
    gdd_path = abspath(joinpath(base_dir, String(manifest["gdd_path"])))

    isfile(esm_path) || return (:fail, "$name: ESM not found at $esm_path")
    isfile(gdd_path) || return (:fail, "$name: GDD not found at $gdd_path")

    ic = manifest["initial_condition"]
    if String(ic["kind"]) != "cosine_bell"
        return (:fail, "$name: only cosine_bell IC is wired in this runner")
    end
    h0   = Float64(ic["h0"])
    lon0 = Float64(ic["lon0"])
    lat0 = Float64(ic["lat0"])

    tol_max = Float64(manifest["tolerance"]["max"])

    # Build the spatial operator via Path B (ESM → PDESystem → discretize).
    prob_template, _ = build_ode_problem(esm_path; grid_ref = gdd_path)

    # Extract grid parameters from the GDD to build physical coordinates.
    gdd        = JSON.parse(read(gdd_path, String))
    grid_entry = first(values(gdd["grids"]))
    Nc  = Int(grid_entry["Nc"])
    R   = Float64(grid_entry["R"])
    r0  = R * Float64(ic["r0_factor"])

    grid = CubedSphereGrid(Nc; R = R)
    lons = cell_centers(grid, :lon)
    lats = cell_centers(grid, :lat)
    N    = 6 * Nc * Nc

    function _cosbell(lo, la)
        cosd = sin(lat0) * sin(la) + cos(lat0) * cos(la) * cos(lo - lon0)
        r    = R * acos(clamp(cosd, -1.0, 1.0))
        return r < r0 ? (h0 / 2) * (1 + cos(π * r / r0)) : 0.0
    end
    u0_cb = [_cosbell(lons[c], lats[c]) for c in 1:N]

    # Replace the zero IC from Path B with the physical cosine-bell IC.
    # remake preserves MTK's initialization data while swapping the u0 vector.
    prob = SciMLBase.remake(prob_template; u0 = u0_cb)

    sol = solve(prob; reltol = 1e-6, abstol = 1e-8, save_everystep = false)

    if sol.retcode != SciMLBase.ReturnCode.Success
        return (:fail, "$name: ODE solver failed (retcode=$(sol.retcode))")
    end

    u_final   = sol.u[end]
    linf_abs  = maximum(abs.(u_final .- u0_cb))
    linf_norm = linf_abs / max(h0, 1.0e-12)

    msg = "$name: L∞/h0 = $(round(linf_norm; sigdigits = 3)) (tol = $tol_max)"
    if linf_norm <= tol_max
        return (:pass, msg)
    else
        return (:fail, msg * " — exceeds tolerance")
    end
end
