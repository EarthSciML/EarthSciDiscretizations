#!/usr/bin/env julia
#=
Regenerate golden d2u/dx2-output arrays for rect_1d_burgers_diffusion_periodic.

Canonical pipeline: load rule JSON → read `replacement` field directly
(rule carries authored replacement; esd-t4h) → evaluate the replacement AST
at each cell via the ESD/ESS passthrough.

Scalar sub-expressions (coefficient terms without `index` ops) are evaluated
through `EarthSciDiscretizations.eval_coeff`, a thin passthrough to
`EarthSciSerialization.evaluate_expr`. Array-indexing ops (`index`) are resolved
against the field vector with periodic wrap. No math is implemented outside the
rule AST — the formula lives entirely in
`discretizations/finite_difference/centered_2nd_deriv_uniform.json`.

Run from the repo root:
    julia --project=. tests/conformance/discretization/rect_1d_burgers_diffusion_periodic/regenerate_golden.jl

SHA header in the golden: MOL #531 at 35cc9143dc553ac7d3619738bd77b250c1ed162f
(conditions taken from MOL #531 Test 00 advection IC — same Gaussian used for
cross-operator consistency with the advection goldens in esd-11y).
=#

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "..")); io = devnull)
    Pkg.add("JSON"; io = devnull)
end

using EarthSciDiscretizations: eval_coeff
using JSON

const HERE = @__DIR__
const ROOT = abspath(joinpath(HERE, "..", "..", "..", ".."))
const FX_PATH = joinpath(HERE, "fixtures.json")
const GLD_DIR = joinpath(HERE, "golden")

const MOL531_SHA = "35cc9143dc553ac7d3619738bd77b250c1ed162f"

# ---------------------------------------------------------------------------
# Recursive replacement-AST evaluator (1D periodic, same logic as advection).
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
    op = String(node["op"])
    args = node["args"]

    if op == "index"
        function _resolve_idx(ie)
            ie isa AbstractString && String(ie) == "\$x" && return cell_idx
            ie isa AbstractDict || error("unsupported index expression: $ie")
            ia = ie["args"]
            String(ia[1]) == "\$x" || error("index lhs not '\$x': $(ia[1])")
            op2 = String(ie["op"])
            op2 == "+" && return cell_idx + Int(ia[2])
            op2 == "-" && return cell_idx - Int(ia[2])
            error("unsupported index op '$op2'")
        end
        j = _resolve_idx(args[2])
        return u[mod(j - 1, N) + 1]
    end

    !_has_index(node) && return eval_coeff(node, bindings)

    ev(a) = _eval_replacement(a, u, cell_idx, bindings, N)
    op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
    op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
    op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
    op == "/" && return ev(args[1]) / ev(args[2])
    error("_eval_replacement: unsupported op '$op'")
end

# ---------------------------------------------------------------------------
# Field construction
# ---------------------------------------------------------------------------

function _gaussian_field(
        n::Int, x_start::Float64, x_end::Float64,
        center::Float64, sigma::Float64, scale::Float64
    )
    dx = (x_end - x_start) / n
    amp = scale / (sigma * sqrt(2π))
    return Float64[
        amp * exp(-(x_start + (i - 0.5) * dx - center)^2 / (2 * sigma^2))
            for i in 1:n
    ]
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    spec = JSON.parsefile(FX_PATH)
    rule_path = joinpath(ROOT, spec["rule_path"])
    raw_rule = JSON.parsefile(rule_path)
    rule_name = spec["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    repl_raw = rule_body["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
        repl_raw["expr"] : repl_raw

    isdir(GLD_DIR) || mkpath(GLD_DIR)

    for fx in spec["fixtures"]
        g = fx["grid"]
        ic = fx["initial_condition"]
        N = Int(g["n_cells"])
        dx = (Float64(g["x_end"]) - Float64(g["x_start"])) / N
        bc = String(g["bc"])
        @assert bc == "periodic" "Only periodic BC supported in this script"

        u = _gaussian_field(
            N, Float64(g["x_start"]), Float64(g["x_end"]),
            Float64(ic["center"]), Float64(ic["sigma"]),
            Float64(ic["scale_factor"])
        )

        bindings = Dict{String, Float64}("dx" => dx)
        d2u_dx2 = [_eval_replacement(replacement, u, i, bindings, N) for i in 1:N]

        out = Dict{String, Any}(
            "_mol531_sha" => MOL531_SHA,
            "_captured_by" => "esd-i73",
            "version" => spec["version"],
            "rule" => rule_name,
            "fixture" => fx["name"],
            "grid" => Dict{String, Any}("n_cells" => N, "dx" => dx, "bc" => bc),
            "field_u" => u,
            "d2u_dx2" => d2u_dx2,
        )

        path = joinpath(GLD_DIR, "$(fx["name"]).json")
        open(path, "w") do io
            JSON.print(io, out, 2)
            write(io, "\n")
        end
        println("wrote $path")
    end
    return
end

main()
