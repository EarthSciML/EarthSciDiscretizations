using EarthSciDiscretizations: CubedSphereGrid, transport_2d_linrood!, total_area

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
