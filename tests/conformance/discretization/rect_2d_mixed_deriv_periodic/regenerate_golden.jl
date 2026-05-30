#!/usr/bin/env julia
#=
Regenerate golden mixed-derivative output arrays for rect_2d_mixed_deriv_periodic.

Canonical pipeline: load rule JSON → evaluate replacement AST at each cell
via the ESD/ESS passthrough.

The rule `mixed_deriv_2nd_uniform` uses a direct `replacement` (no stencil),
so `lower_stencil_to_replacement` is a no-op (idempotent). The replacement
expression is the 4-point stencil:
  (u[i+1,j+1] - u[i+1,j-1] - u[i-1,j+1] + u[i-1,j-1]) / (4*dx*dy)

Scalar sub-expressions without `index` ops are evaluated through
`EarthSciDiscretizations.eval_coeff`, a thin passthrough to
`EarthSciSerialization.evaluate_expr`. Array-indexing ops (`index`) are
resolved against the 2D field matrix with periodic wrap in both dimensions.
No math is implemented outside the rule AST.

Run from the repo root:
    julia --project=. tests/conformance/discretization/rect_2d_mixed_deriv_periodic/regenerate_golden.jl

Field: u(x,y) = sin(x)*sin(y) on [0, 2π]×[0, 2π]. Exact Dxy = cos(x)*cos(y).
SHA header in the golden: MOL #531 at 35cc9143dc953ac7d3619738bd77b250c1ed162f
=#

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "..")); io = devnull)
    Pkg.add("JSON"; io = devnull)
end

using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
using JSON

const HERE    = @__DIR__
const ROOT    = abspath(joinpath(HERE, "..", "..", "..", ".."))
const FX_PATH = joinpath(HERE, "fixtures.json")
const GLD_DIR = joinpath(HERE, "golden")

const MOL531_SHA = "35cc9143dc953ac7d3619738bd77b250c1ed162f"

function _has_index(node)::Bool
    node isa AbstractDict || return false
    node["op"] == "index" && return true
    return any(_has_index(a) for a in get(node, "args", Any[]))
end

function _resolve_axis_idx(ie, xi::Int, yi::Int)::Int
    if ie isa AbstractString
        s = String(ie)
        s == "\$x" && return xi
        s == "\$y" && return yi
        error("_resolve_axis_idx: unknown axis variable '$s'")
    end
    ie isa AbstractDict || error("_resolve_axis_idx: unexpected type $(typeof(ie))")
    ia  = ie["args"]
    op2 = String(ie["op"])
    lhs = String(ia[1])
    off = Int(ia[2])
    base = lhs == "\$x" ? xi : (lhs == "\$y" ? yi : error("unknown lhs '$lhs'"))
    op2 == "+" && return base + off
    op2 == "-" && return base - off
    error("_resolve_axis_idx: unsupported op '$op2'")
end

function _eval_replacement_2d(
        node,
        u::Matrix{Float64},
        xi::Int, yi::Int,
        bindings::Dict{String, Float64},
        Nx::Int, Ny::Int,
)::Float64
    node isa Number && return Float64(node)
    if node isa AbstractString
        s = String(node)
        s == "\$x" && return Float64(xi)
        s == "\$y" && return Float64(yi)
        haskey(bindings, s) && return bindings[s]
        error("_eval_replacement_2d: unresolved variable '$s'")
    end
    node isa AbstractDict || error("unexpected node type $(typeof(node))")
    op   = String(node["op"])
    args = node["args"]
    if op == "index"
        length(args) == 3 || error("index expects 3 args, got $(length(args))")
        jx = _resolve_axis_idx(args[2], xi, yi)
        jy = _resolve_axis_idx(args[3], xi, yi)
        return u[mod(jx - 1, Nx) + 1, mod(jy - 1, Ny) + 1]
    end
    !_has_index(node) && return eval_coeff(node, bindings)
    ev(a) = _eval_replacement_2d(a, u, xi, yi, bindings, Nx, Ny)
    op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
    op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
    op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
    op == "/" && return ev(args[1]) / ev(args[2])
    error("unsupported op '$op'")
end

rule_path = joinpath(ROOT, "discretizations", "finite_difference",
                     "mixed_deriv_2nd_uniform.json")
raw_rule  = JSON.parsefile(rule_path)
rule_body = raw_rule["discretizations"]["mixed_deriv_2nd_uniform"]
rule_lowered = lower_stencil_to_replacement(rule_body)
repl_raw  = rule_lowered["replacement"]
replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
              repl_raw["expr"] : repl_raw

fxs = JSON.parsefile(FX_PATH)
mkpath(GLD_DIR)

for fx in fxs["fixtures"]
    g   = fx["grid"]
    Nx  = Int(g["n_cells_x"])
    Ny  = Int(g["n_cells_y"])
    x0  = Float64(g["x_start"])
    x1  = Float64(g["x_end"])
    y0  = Float64(g["y_start"])
    y1  = Float64(g["y_end"])
    dx  = (x1 - x0) / Nx
    dy  = (y1 - y0) / Ny
    xs  = [x0 + (i - 0.5) * dx for i in 1:Nx]
    ys  = [y0 + (j - 0.5) * dy for j in 1:Ny]

    # u(x,y) = sin(x) * sin(y)
    u = [sin(xs[i]) * sin(ys[j]) for i in 1:Nx, j in 1:Ny]

    bindings = Dict{String, Float64}("dx" => dx, "dy" => dy)

    dxy_u = Matrix{Float64}(undef, Nx, Ny)
    for iy in 1:Ny, ix in 1:Nx
        dxy_u[ix, iy] = _eval_replacement_2d(replacement, u, ix, iy, bindings, Nx, Ny)
    end

    field_rows = [[dxy_u[i, j] for j in 1:Ny] for i in 1:Nx]

    out = Dict(
        "_captured_by" => "esd-wdv",
        "_mol531_sha"  => MOL531_SHA,
        "rule"         => "mixed_deriv_2nd_uniform",
        "fixture"      => fx["name"],
        "grid"         => Dict(
            "bc"        => "periodic",
            "dx"        => dx,
            "dy"        => dy,
            "n_cells_x" => Nx,
            "n_cells_y" => Ny,
        ),
        "dxy_u"        => field_rows,
    )

    out_path = joinpath(GLD_DIR, fx["name"] * ".json")
    open(out_path, "w") do io
        JSON.print(io, out)
    end
    println("Generated: ", relpath(out_path, ROOT))
end
