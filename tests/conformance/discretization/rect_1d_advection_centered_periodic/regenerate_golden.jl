#!/usr/bin/env julia
#=
Regenerate golden grad-output arrays for rect_1d_advection_centered_periodic.

Canonical pipeline: load rule JSON → lower_stencil_to_replacement (idempotent
for centered_2nd_uniform, which already carries a `replacement` field) → evaluate
the replacement AST at each cell via the ESD/ESS passthrough.

Scalar sub-expressions (coefficient terms without `index` ops) are evaluated
through `EarthSciDiscretizations.eval_coeff`, a thin passthrough to
`EarthSciSerialization.evaluate_expr`. Array-indexing ops (`index`) are resolved
against the field vector with periodic wrap. No math is implemented outside the
rule AST — the formula lives entirely in
`discretizations/finite_difference/centered_2nd_uniform.json`.

Run from the repo root:
    julia --project=. tests/conformance/discretization/rect_1d_advection_centered_periodic/regenerate_golden.jl

SHA header in the golden: MOL #531 at 35cc9143dc553ac7d3619738bd77b250c1ed162f
(conditions taken from Test 00 of test/pde_systems/MOL_1D_Linear_Convection.jl).
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

const MOL531_SHA = "35cc9143dc553ac7d3619738bd77b250c1ed162f"

# ---------------------------------------------------------------------------
# Recursive replacement-AST evaluator.
#
# Traverses the replacement Dict produced by `lower_stencil_to_replacement`.
# Pure-scalar sub-trees (no `index` ops) delegate to `eval_coeff` which routes
# through EarthSciSerialization.evaluate_expr — the canonical ESS scalar path.
# Mixed or `index` nodes are recursed; array lookups use 1-based periodic wrap.
# ---------------------------------------------------------------------------

function _has_index(node)::Bool
    node isa AbstractDict || return false
    node["op"] == "index" && return true
    args = get(node, "args", Any[])
    return any(_has_index(a) for a in args)
end

"""
    _eval_replacement(node, u, cell_idx, bindings, N)

Evaluate one node of the replacement AST.
- `u`         : 1-based Float64 vector, length N
- `cell_idx`  : 1-based current cell index (bound to `\$x`)
- `bindings`  : scalar bindings forwarded to `eval_coeff` (e.g. `{"dx" => dx}`)
- `N`         : number of cells (for periodic wrap)
"""
function _eval_replacement(
        node,
        u::Vector{Float64},
        cell_idx::Int,
        bindings::Dict{String, Float64},
        N::Int,
)::Float64
    # ---- literal / scalar variable ----------------------------------------
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

    # ---- index op ----------------------------------------------------------
    if op == "index"
        # args[1] is the operand name ("\$u"); args[2] is the index expression
        function _resolve_idx(ie)
            ie isa AbstractString && String(ie) == "\$x" && return cell_idx
            ie isa AbstractDict || error("_eval_replacement: unsupported index: $ie")
            ia = ie["args"]
            String(ia[1]) == "\$x" ||
                error("_eval_replacement: index lhs is not '\$x': $(ia[1])")
            op2 = String(ie["op"])
            op2 == "+" && return cell_idx + Int(ia[2])
            op2 == "-" && return cell_idx - Int(ia[2])
            error("_eval_replacement: unsupported index op '$op2'")
        end
        j = _resolve_idx(args[2])
        return u[mod(j - 1, N) + 1]   # 1-based periodic wrap
    end

    # ---- pure-scalar sub-tree → eval_coeff (ESS canonical path) ------------
    if !_has_index(node)
        return eval_coeff(node, bindings)
    end

    # ---- mixed ops: recurse ------------------------------------------------
    ev(a) = _eval_replacement(a, u, cell_idx, bindings, N)
    if op == "+"
        length(args) == 1 && return ev(args[1])
        return sum(ev(a) for a in args)
    elseif op == "-"
        length(args) == 1 && return -ev(args[1])
        length(args) == 2 && return ev(args[1]) - ev(args[2])
    elseif op == "*"
        length(args) == 1 && return ev(args[1])
        return prod(ev(a) for a in args)
    elseif op == "/"
        return ev(args[1]) / ev(args[2])
    end
    error("_eval_replacement: unsupported op '$op' in mixed expression")
end

# ---------------------------------------------------------------------------
# Field construction
# ---------------------------------------------------------------------------

function _gaussian_field(n::Int, x_start::Float64, x_end::Float64,
                         center::Float64, sigma::Float64, scale::Float64)
    dx  = (x_end - x_start) / n
    amp = scale / (sigma * sqrt(2π))
    return Float64[amp * exp(-(x_start + (i - 0.5) * dx - center)^2 / (2 * sigma^2))
                   for i in 1:n]
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
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw = rule_lowered["replacement"]
    # centered_2nd_uniform ships its replacement wrapped in {op:"arrayop",...};
    # the arithmetic expression lives in repl_raw["expr"].
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

        u = _gaussian_field(N, Float64(g["x_start"]), Float64(g["x_end"]),
                            Float64(ic["center"]), Float64(ic["sigma"]),
                            Float64(ic["scale_factor"]))

        bindings = Dict{String, Float64}("dx" => dx)
        du_dx = [_eval_replacement(replacement, u, i, bindings, N) for i in 1:N]

        out = Dict{String, Any}(
            "_mol531_sha"  => MOL531_SHA,
            "_captured_by" => "esd-11y",
            "version"      => spec["version"],
            "rule"         => rule_name,
            "fixture"      => fx["name"],
            "grid"         => Dict{String, Any}("n_cells" => N, "dx" => dx, "bc" => bc),
            "field_u"      => u,
            "du_dx"        => du_dx,
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
