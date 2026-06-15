using EarthSciDiscretizations: CubedSphereGrid, total_area, build_ode_problem, cell_centers
import SciMLBase
using OrdinaryDiffEqDefault: solve

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
