#!/usr/bin/env julia
# Regenerate the golden for the point cell-averaging regridder conformance harness.
#
# The golden is the ANALYTIC oracle (no engine): each cell is assigned the integer
# spatial bin (floor(lon/dx), floor(lat/dy)); the contributing stations are those
# sharing the cell's bin; F_tgt = their arithmetic mean, or the configured
# missing_value when the bin is empty. count_denom is the host-precomputed
# contributing-station tally the engine's group-by must reproduce (it is also the
# average's denominator), so it is emitted into the golden as the engine input.
# The reference binding's conformance test (test/test_regridding_conformance.jl)
# feeds station/cell/val/count_denom through the ESS value-invention front-door and
# checks count + F_tgt against this oracle.
#
# Run:  julia tests/conformance/regridding/point/regenerate_golden.jl

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add("JSON"; io = devnull)
end
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

bin(coord, d) = floor(Int, coord / d)

function main()
    spec = JSON.parsefile(FIXTURES)
    dx = Float64(spec["dx"]); dy = Float64(spec["dy"])
    missing_value = Float64(spec["missing_value"])
    slon = Float64.(spec["stations"]["lon"]); slat = Float64.(spec["stations"]["lat"])
    sval = Float64.(spec["stations"]["val"])
    clon = Float64.(spec["cells"]["lon"]); clat = Float64.(spec["cells"]["lat"])

    nstation = length(slon); ncell = length(clon)
    sbin = [(bin(slon[i], dx), bin(slat[i], dy)) for i in 1:nstation]
    cbin = [(bin(clon[j], dx), bin(clat[j], dy)) for j in 1:ncell]

    count_denom = zeros(Float64, ncell)
    F_tgt = zeros(Float64, ncell)
    for j in 1:ncell
        contrib = [i for i in 1:nstation if sbin[i] == cbin[j]]
        count_denom[j] = length(contrib)
        F_tgt[j] = isempty(contrib) ? missing_value : sum(sval[contrib]) / length(contrib)
    end

    out = Dict{String, Any}(
        "model" => "point",
        "reference_binding" => "julia",
        "missing_value" => missing_value,
        "count_denom" => count_denom,   # engine input (host tally) AND expected count
        "count" => count_denom,         # expected engine group-by tally
        "F_tgt" => F_tgt,               # expected per-cell mean / missing_value
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN  (count_denom=$count_denom, F_tgt=$F_tgt)")
    return
end

main()
