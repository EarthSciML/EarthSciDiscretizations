#!/usr/bin/env julia
# Regenerate the golden for the lev=min 3D->surface reduction conformance harness.
#
# The golden is the ANALYTIC oracle (no engine): the surface field is the F_3d
# slice at the level whose vertical coordinate lev_coord is the unique minimum
# (an argmin gather). F_surf[x][y] = F_3d[x][y][argmin_lev]. This is the
# closed-form reduction the ESS engine's declarative
# `lev_min_surface_reduce.esm` must reproduce (sum_product over lev with an
# inline min_sum equality indicator). The reference binding's conformance test
# (test/test_reduction_conformance.jl) feeds F_3d + lev_coord through
# build_evaluator and checks F_surf against this oracle.
#
# Run:  julia tests/conformance/reduction/lev_min/regenerate_golden.jl

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add("JSON"; io = devnull)
end
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

function main()
    spec = JSON.parsefile(FIXTURES)
    lev_coord = Float64.(spec["lev_coord"])
    F_3d = spec["F_3d"]                       # nested [x][y][lev]
    nx = length(F_3d); ny = length(F_3d[1]); nlev = length(F_3d[1][1])

    # Unique minimum of the vertical coordinate (the value:min level).
    lmin = argmin(lev_coord)
    @assert count(==(lev_coord[lmin]), lev_coord) == 1 "lev_coord must have a UNIQUE minimum"

    # Surface = the F_3d slice at the argmin level (F_surf[x][y]).
    F_surf = [[Float64(F_3d[a][b][lmin]) for b in 1:ny] for a in 1:nx]

    out = Dict{String, Any}(
        "model" => "lev_min",
        "reference_binding" => "julia",
        "argmin_level" => lmin,        # 1-based level index of the minimum coordinate
        "F_surf" => F_surf,            # expected surfaced field, F_surf[x][y]
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN  (argmin_level=$lmin, F_surf=$F_surf)")
    return
end

main()
