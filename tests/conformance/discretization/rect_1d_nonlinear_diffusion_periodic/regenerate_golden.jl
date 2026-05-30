#!/usr/bin/env julia
#=
Regenerate golden nonlinear-Laplacian-output arrays for
rect_1d_nonlinear_diffusion_periodic.

Canonical pipeline: load rule JSON (nonlinear_laplacian_uniform already ships in
`replacement`/arrayop form — no stencil lowering needed) → evaluate the replacement
AST at each cell via the ESD/ESS passthrough.

The rule computes d_x(f * d_x(u)) discretized as:
  [(f[x]+f[x+1])/2*(u[x+1]-u[x]) - (f[x-1]+f[x])/2*(u[x]-u[x-1])] / dx²

Both fields $u and $f are resolved against separate 1D vectors with periodic wrap.
Scalar sub-expressions delegate to `eval_coeff` (ESS canonical path).

Run from the repo root:
    julia --project=. tests/conformance/discretization/rect_1d_nonlinear_diffusion_periodic/regenerate_golden.jl

SHA header in the golden: MOL #531 at 35cc9143dc553ac7d3619738bd77b250c1ed162f
(conditions derived from MOL #531 NonLinear_Diffusion Test 00, adapted for periodic
BCs using smooth sinusoidal fields).
=#

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "..")); io = devnull)
    Pkg.add("JSON"; io = devnull)
end

using EarthSciDiscretizations: eval_coeff
using JSON

const HERE    = @__DIR__
const ROOT    = abspath(joinpath(HERE, "..", "..", "..", ".."))
const FX_PATH = joinpath(HERE, "fixtures.json")
const GLD_DIR = joinpath(HERE, "golden")

const MOL531_SHA = "35cc9143dc553ac7d3619738bd77b250c1ed162f"

# ---------------------------------------------------------------------------
# Recursive replacement-AST evaluator for two-field 1D periodic case.
#
# The nonlinear_laplacian_uniform replacement references two fields: $u and $f.
# Index ops carry the field name as args[1] ("$u" or "$f"), and a 1D axis
# expression as args[2] (either "$x" or an offset dict).
# ---------------------------------------------------------------------------

function _has_index(node)::Bool
    node isa AbstractDict || return false
    node["op"] == "index" && return true
    args = get(node, "args", Any[])
    return any(_has_index(a) for a in args)
end

function _eval_replacement(
        node,
        u::Vector{Float64},
        f::Vector{Float64},
        cell_idx::Int,
        bindings::Dict{String, Float64},
        N::Int,
)::Float64
    node isa Number && return Float64(node)
    if node isa AbstractString
        s = String(node)
        s == "\$x" && return Float64(cell_idx)
        haskey(bindings, s) && return bindings[s]
        error("_eval_replacement: unresolved variable '$s'")
    end

    node isa AbstractDict || error("_eval_replacement: unexpected node type $(typeof(node))")
    op   = String(node["op"])
    args = node["args"]

    if op == "index"
        # args[1] is the field name ("\$u" or "\$f"), args[2] is the index expression
        field_name = String(args[1])
        field = field_name == "\$u" ? u :
                field_name == "\$f" ? f :
                error("_eval_replacement: unknown field '$field_name'")
        function _resolve_idx(ie)
            ie isa AbstractString && String(ie) == "\$x" && return cell_idx
            ie isa AbstractDict || error("unsupported index expression: $ie")
            ia  = ie["args"]
            String(ia[1]) == "\$x" || error("index lhs not '\$x': $(ia[1])")
            op2 = String(ie["op"])
            op2 == "+" && return cell_idx + Int(ia[2])
            op2 == "-" && return cell_idx - Int(ia[2])
            error("unsupported index op '$op2'")
        end
        j = _resolve_idx(args[2])
        return field[mod(j - 1, N) + 1]
    end

    !_has_index(node) && return eval_coeff(node, bindings)

    ev(a) = _eval_replacement(a, u, f, cell_idx, bindings, N)
    op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
    op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
    op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
    op == "/" && return ev(args[1]) / ev(args[2])
    error("_eval_replacement: unsupported op '$op'")
end

# ---------------------------------------------------------------------------
# Field construction
# ---------------------------------------------------------------------------

function _sinusoidal_plus_const_field(n::Int, x_start::Float64, x_end::Float64,
                                      amplitude::Float64, constant::Float64)
    dx = (x_end - x_start) / n
    return Float64[constant + amplitude * sin(π * (x_start + (i - 0.5) * dx))
                   for i in 1:n]
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    spec = JSON.parsefile(FX_PATH)
    rule_path = joinpath(ROOT, spec["rule_path"])
    raw_rule  = JSON.parsefile(rule_path)
    rule_name = spec["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    # nonlinear_laplacian_uniform already ships with a `replacement` field
    repl_raw  = rule_body["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    isdir(GLD_DIR) || mkpath(GLD_DIR)

    for fx in spec["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        N   = Int(g["n_cells"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / N
        bc  = String(g["bc"])
        @assert bc == "periodic" "Only periodic BC supported in this script"

        u_ic = ic["u"]
        u = _sinusoidal_plus_const_field(
            N, Float64(g["x_start"]), Float64(g["x_end"]),
            Float64(u_ic["amplitude"]), Float64(u_ic["constant"]),
        )
        f = 1.0 ./ u   # f(x) = 1/u(x), matches fixtures.json formula

        bindings = Dict{String, Float64}("dx" => dx)
        nl_lap_u = [_eval_replacement(replacement, u, f, i, bindings, N) for i in 1:N]

        out = Dict{String, Any}(
            "_mol531_sha"        => MOL531_SHA,
            "_captured_by"       => "esd-i73",
            "version"            => spec["version"],
            "rule"               => rule_name,
            "fixture"            => fx["name"],
            "grid"               => Dict{String, Any}("n_cells" => N, "dx" => dx, "bc" => bc),
            "field_u"            => u,
            "field_f"            => f,
            "nonlinear_laplacian_u" => nl_lap_u,
        )

        path = joinpath(GLD_DIR, "$(fx["name"]).json")
        open(path, "w") do io
            JSON.print(io, out, 2)
            write(io, "\n")
        end
        println("wrote $path")
    end
end

main()
