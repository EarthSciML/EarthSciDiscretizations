#!/usr/bin/env julia
# Regenerate the golden for the staggered B-spline regridder conformance harness.
#
# The golden is an INDEPENDENT closed-form ORACLE (no engine): the host-derived
# kernel inputs (F_src = the test polynomial sampled at the uniform source nodes;
# base = the per-target stencil-start index; s = the per-target local coordinate)
# plus the expected F_tgt = the SAME polynomial evaluated directly at the target
# locations. Polynomial reproduction is the B-spline invariant (RFC §8): a degree-d
# kernel fed a degree-≤d field must return the field's value at the target. The
# reference binding's conformance test (test/test_regridding_conformance.jl) feeds
# F_src/base/s through the ESS engine and checks the engine reproduces this oracle.
#
# Run:  julia tests/conformance/regridding/bspline/regenerate_golden.jl

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add("JSON"; io = devnull)
end
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

poly1d(coeffs, x) = sum(Float64(coeffs[k]) * x^(k - 1) for k in 1:length(coeffs))
poly2d(c, x, y) = Float64(c["c00"]) + Float64(c["cx"]) * x + Float64(c["cy"]) * y +
    Float64(c["cxy"]) * x * y

base_1d(rule, x) = rule == "floor(x) + 1" ? (floor(Int, x) + 1) :
    rule == "floor(x)" ? floor(Int, x) :
    error("unknown base_rule $rule")
frac(x) = x - floor(x)

function build_model(m)
    name = m["name"]
    out = Dict{String, Any}("name" => name)
    if m["ndim"] == 1
        nodes = Int(m["nodes"])
        coeffs = m["polynomial"]["coeffs"]
        xs = Float64.(0:(nodes - 1))
        out["F_src"] = [poly1d(coeffs, x) for x in xs]
        tx = Float64.(m["targets_x"])
        out["base"] = [base_1d(m["base_rule"], x) for x in tx]
        out["s"] = [frac(x) for x in tx]
        out["F_tgt"] = [poly1d(coeffs, x) for x in tx]   # ORACLE: poly at targets
    elseif m["ndim"] == 2
        nx = Int(m["nodes_x"]); ny = Int(m["nodes_y"])
        c = m["polynomial"]["coeffs"]
        xs = Float64.(0:(nx - 1)); ys = Float64.(0:(ny - 1))
        # F_src[i][j] = P(x_i, y_j) — reconstructed as a Matrix [i, j] by the consumer.
        out["F_src"] = [[poly2d(c, x, y) for y in ys] for x in xs]
        tgt = m["targets_xy"]
        out["base_x"] = [floor(Int, Float64(t[1])) + 1 for t in tgt]
        out["base_y"] = [floor(Int, Float64(t[2])) + 1 for t in tgt]
        out["s_x"] = [frac(Float64(t[1])) for t in tgt]
        out["s_y"] = [frac(Float64(t[2])) for t in tgt]
        out["F_tgt"] = [poly2d(c, Float64(t[1]), Float64(t[2])) for t in tgt]
    else
        error("unsupported ndim $(m["ndim"])")
    end
    return out
end

function main()
    spec = JSON.parsefile(FIXTURES)
    models = [build_model(m) for m in spec["models"]]
    out = Dict{String, Any}(
        "model" => "bspline",
        "reference_binding" => "julia",
        "base_indexing" => spec["base_indexing"],
        "models" => models,
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN")
    return
end

main()
