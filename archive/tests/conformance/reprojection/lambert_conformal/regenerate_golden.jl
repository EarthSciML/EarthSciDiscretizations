#!/usr/bin/env julia
# Regenerate the golden for the lambert_conformal reprojection conformance harness.
#
# golden.forward is an INDEPENDENT closed-form ORACLE: the spherical Lambert
# Conformal Conic forward map (Snyder 1987 eqs. 15-1..15-4, identical to PROJ
# `+proj=lcc +R=...` for the sphere), implemented here directly — NOT by evaluating
# the .esm AST. The reference binding's conformance test
# (test/test_reprojection_conformance.jl) then drives reprojection/lambert_conformal.esm
# through the ESS engine and checks it reproduces this oracle to the declared
# tolerance, so an independent implementation cross-checks the declarative document.
#
# As a guard, this script asserts the oracle reproduces the proj4-validated
# reference points (lifted from reprojection/lambert_conformal_conformance_test.jl,
# bead esd-47z.2) before writing — if the closed form is mistranscribed the
# regenerate fails loudly rather than freezing a wrong golden.
#
# Run:  julia tests/conformance/reprojection/lambert_conformal/regenerate_golden.jl

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add("JSON"; io = devnull)
end
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

deg2rad_(d) = d * (pi / 180.0)

# Spherical LCC forward, Snyder 1987 (two standard parallels, lat_1 != lat_2).
# inputs in degrees, x/y in metres.
function lcc_forward(lon, lat; lat_1, lat_2, lat_0, lon_0, R)
    phi1 = deg2rad_(lat_1); phi2 = deg2rad_(lat_2)
    phi0 = deg2rad_(lat_0); phi = deg2rad_(lat)
    lam0 = deg2rad_(lon_0); lam = deg2rad_(lon)
    n = log(cos(phi1) / cos(phi2)) /
        log(tan(pi / 4 + phi2 / 2) / tan(pi / 4 + phi1 / 2))
    F = cos(phi1) * tan(pi / 4 + phi1 / 2)^n / n
    rho = R * F / tan(pi / 4 + phi / 2)^n
    rho0 = R * F / tan(pi / 4 + phi0 / 2)^n
    x = rho * sin(n * (lam - lam0))
    y = rho0 - rho * cos(n * (lam - lam0))
    return (x, y)
end

# proj4-validated reference forward (esd-47z.2): {param_set => [(lon,lat,x,y)...]}.
const _GOLD_REF = Dict(
    "WRF" => [
        (-97.0, 39.0, 0.0, 0.43226828519254923),
        (-120.0, 35.0, -2028208.5469169612, -140947.2885132283),
        (-75.0, 40.0, 1795192.9893785133, 356152.97854254674),
        (-90.0, 45.0, 530758.0790195653, 668954.7596800169),
        (-100.0, 25.0, -309853.27839640283, -1541532.7554675443),
        (-80.0, 48.0, 1213070.0930733178, 1097145.5185974706),
        (-110.0, 31.0, -1228247.9439207606, -773888.551179463),
        (-104.0, 42.5, -554207.2432213794, 401411.8990695188),
    ],
    "NEI2016" => [
        (-97.0, 39.0, 0.0, -110589.55965320487),
        (-120.0, 35.0, -2066463.0536515424, -290403.0815928243),
        (-75.0, 40.0, 1845776.353572973, 224515.92418459523),
        (-90.0, 45.0, 549842.4483977047, 575299.4792623082),
        (-100.0, 25.0, -309377.0930418335, -1668877.3176074345),
        (-80.0, 48.0, 1266623.2685378059, 1007635.9929630421),
        (-110.0, 31.0, -1239963.1270495378, -909332.8255445026),
        (-104.0, 42.5, -571190.8151147809, 298694.874228091),
    ],
)

isclose(a, b; rtol = 1.0e-7, atol = 1.0e-4) = abs(a - b) <= atol + rtol * abs(b)

function main()
    spec = JSON.parsefile(FIXTURES)
    coords = [(Float64(c[1]), Float64(c[2])) for c in spec["coords"]]
    out_sets = Any[]
    for ps in spec["param_sets"]
        name = ps["name"]
        p = ps["params"]
        kw = (
            lat_1 = Float64(p["lat_1"]), lat_2 = Float64(p["lat_2"]),
            lat_0 = Float64(p["lat_0"]), lon_0 = Float64(p["lon_0"]), R = Float64(p["R"]),
        )
        fwd = Vector{Vector{Float64}}()
        for (lon, lat) in coords
            x, y = lcc_forward(lon, lat; kw...)
            push!(fwd, [x, y])
        end
        # GUARD: the oracle must reproduce the proj4-validated reference.
        ref = _GOLD_REF[name]
        @assert length(ref) == length(coords) "ref/coords length mismatch for $name"
        for (k, (rlon, rlat, rx, ry)) in enumerate(ref)
            @assert coords[k][1] == rlon && coords[k][2] == rlat "coord mismatch $name k=$k"
            isclose(fwd[k][1], rx) || error("$name pt $k: x oracle $(fwd[k][1]) != proj4 $rx")
            isclose(fwd[k][2], ry) || error("$name pt $k: y oracle $(fwd[k][2]) != proj4 $ry")
        end
        push!(
            out_sets, Dict{String, Any}(
                "name" => name,
                "params" => p,
                "forward" => fwd,
                "roundtrip" => [collect(c) for c in coords],   # inverse-of-forward = identity
            )
        )
        println("  $name: oracle reproduces proj4 reference (8 pts) ✓")
    end
    out = Dict{String, Any}(
        "model" => "lambert_conformal",
        "reference_binding" => "julia",
        "coords" => [collect(c) for c in coords],
        "param_sets" => out_sets,
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN")
    return
end

main()
