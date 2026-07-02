#!/usr/bin/env julia
# Regenerate the golden for the longlat reprojection conformance harness.
#
# longlat is the geographic WGS84 identity projection, so the golden is an
# INDEPENDENT CLOSED-FORM ORACLE (no engine, no binding): forward x = lon - lon_0,
# forward y = lat; inverse-of-forward = identity. The reference binding's
# conformance test (test/test_reprojection_conformance.jl) then drives
# reprojection/longlat.esm through the ESS engine and checks the engine reproduces
# this oracle to the declared tolerance — proving the declarative document is green
# through the ESS engine and the golden is consumable by every binding.
#
# Run from anywhere:  julia tests/conformance/reprojection/longlat/regenerate_golden.jl
# Only JSON is needed; activate a throwaway env so the repo deps are untouched.

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add("JSON"; io = devnull)
end
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

# FORWARD map:  x = lon - lon_0,  y = lat.   INVERSE-of-FORWARD = identity.
forward(lon, lat, lon_0) = (lon - lon_0, lat)
roundtrip(lon, lat, lon_0) = let (x, y) = forward(lon, lat, lon_0)
    (x + lon_0, y)   # geo_lon = proj_x + lon_0; geo_lat = proj_y
end

function main()
    spec = JSON.parsefile(FIXTURES)
    coords = [(Float64(c[1]), Float64(c[2])) for c in spec["coords"]]
    cases = Any[]
    for case in spec["cases"]
        lon_0 = Float64(case["lon_0"])
        fwd = [collect(forward(lon, lat, lon_0)) for (lon, lat) in coords]
        rt = [collect(roundtrip(lon, lat, lon_0)) for (lon, lat) in coords]
        push!(
            cases, Dict{String, Any}(
                "name" => case["name"],
                "lon_0" => lon_0,
                "forward" => fwd,
                "roundtrip" => rt,
            )
        )
    end
    out = Dict{String, Any}(
        "model" => "longlat",
        "reference_binding" => "julia",
        "coords" => [collect(c) for c in coords],
        "cases" => cases,
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN")
    return
end

main()
