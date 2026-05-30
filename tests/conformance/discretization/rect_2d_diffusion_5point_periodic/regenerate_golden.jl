#!/usr/bin/env julia
#=
Regenerate golden Laplacian-output arrays for rect_2d_diffusion_5point_periodic.

Canonical pipeline: load rule JSON → lower_stencil_to_replacement (arakawa stencil
with 6 entries over $x and $y axes) → evaluate the replacement AST at each cell
via the ESD/ESS passthrough.

Scalar sub-expressions (coefficient terms without `index` ops) are evaluated
through `EarthSciDiscretizations.eval_coeff`, a thin passthrough to
`EarthSciSerialization.evaluate_expr`. Array-indexing ops (`index`) are resolved
against the 2D field matrix with periodic wrap in both dimensions. No math is
implemented outside the rule AST — the formula lives entirely in
`discretizations/finite_difference/laplacian_2nd_uniform_cartesian.json`.

Run from the repo root:
    julia --project=. tests/conformance/discretization/rect_2d_diffusion_5point_periodic/regenerate_golden.jl

SHA header in the golden: MOL #531 at 35cc9143dc553ac7d3619738bd77b250c1ed162f
(conditions derived from MOL #531 2D diffusion test, adapted for periodic BCs).
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
# Recursive replacement-AST evaluator for 2D periodic fields.
#
# The arakawa-lowered stencil produces index ops with 3 args:
#   {"op": "index", "args": ["$u", x_arg, y_arg]}
# where x_arg and y_arg are either a bare axis variable ("$x"/"$y") or an
# offset expression ({"op": "+", "args": ["$x", k]}).
#
# Pure-scalar sub-trees (no `index` ops) delegate to `eval_coeff` which routes
# through EarthSciSerialization.evaluate_expr — the canonical ESS scalar path.
# ---------------------------------------------------------------------------

function _has_index(node)::Bool
    node isa AbstractDict || return false
    node["op"] == "index" && return true
    args = get(node, "args", Any[])
    return any(_has_index(a) for a in args)
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

"""
    _eval_replacement_2d(node, u, xi, yi, bindings, Nx, Ny)

Evaluate one node of the 2D replacement AST.
- `u`       : Nx×Ny Float64 matrix (1-based)
- `xi`/`yi` : 1-based current cell indices (bound to \$x and \$y)
- `bindings`: scalar bindings forwarded to `eval_coeff` (e.g. dx, dy)
- `Nx`/`Ny` : grid dimensions for periodic wrap
"""
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
        # args: [field_name, x_arg, y_arg]
        length(args) == 3 || error("index op expects 3 args, got $(length(args))")
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
    error("_eval_replacement_2d: unsupported op '$op'")
end

# ---------------------------------------------------------------------------
# Field construction
# ---------------------------------------------------------------------------

function _gaussian_2d_field(nx::Int, ny::Int,
                             x_start::Float64, x_end::Float64,
                             y_start::Float64, y_end::Float64,
                             cx::Float64, cy::Float64,
                             sigma::Float64, scale::Float64)
    dx  = (x_end - x_start) / nx
    dy  = (y_end - y_start) / ny
    amp = scale / (2π * sigma^2)
    u   = Matrix{Float64}(undef, nx, ny)
    for ix in 1:nx
        x = x_start + (ix - 0.5) * dx
        for iy in 1:ny
            y = y_start + (iy - 0.5) * dy
            u[ix, iy] = amp * exp(-((x - cx)^2 + (y - cy)^2) / (2 * sigma^2))
        end
    end
    return u
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
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw  = rule_lowered["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    isdir(GLD_DIR) || mkpath(GLD_DIR)

    for fx in spec["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        Nx  = Int(g["n_cells_x"])
        Ny  = Int(g["n_cells_y"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / Nx
        dy  = (Float64(g["y_end"]) - Float64(g["y_start"])) / Ny
        bc  = String(g["bc"])
        @assert bc == "periodic" "Only periodic BC supported in this script"

        u = _gaussian_2d_field(
            Nx, Ny,
            Float64(g["x_start"]), Float64(g["x_end"]),
            Float64(g["y_start"]), Float64(g["y_end"]),
            Float64(ic["center_x"]), Float64(ic["center_y"]),
            Float64(ic["sigma"]), Float64(ic["scale_factor"]),
        )

        bindings = Dict{String, Float64}("dx" => dx, "dy" => dy)

        lap_u = Matrix{Float64}(undef, Nx, Ny)
        for iy in 1:Ny, ix in 1:Nx
            lap_u[ix, iy] = _eval_replacement_2d(replacement, u, ix, iy, bindings, Nx, Ny)
        end

        # Store as array-of-rows: laplacian_u[ix][iy]
        out = Dict{String, Any}(
            "_mol531_sha"  => MOL531_SHA,
            "_captured_by" => "esd-i73",
            "version"      => spec["version"],
            "rule"         => rule_name,
            "fixture"      => fx["name"],
            "grid"         => Dict{String, Any}(
                "n_cells_x" => Nx, "n_cells_y" => Ny,
                "dx" => dx, "dy" => dy, "bc" => bc,
            ),
            "field_u"      => [u[ix, :] for ix in 1:Nx],
            "laplacian_u"  => [lap_u[ix, :] for ix in 1:Nx],
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
