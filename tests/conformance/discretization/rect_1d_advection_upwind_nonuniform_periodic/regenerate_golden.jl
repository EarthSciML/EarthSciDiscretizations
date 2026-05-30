#!/usr/bin/env julia
#=
Regenerate golden values for the rect_1d_advection_upwind_nonuniform_periodic
conformance harness.

Drives the canonical ESD pipeline:
  1. Load upwind_1st_nonuniform rule from discretizations/
  2. Lower stencil form to replacement via lower_stencil_to_replacement
  3. Evaluate per-cell using the lowered AST with per-cell dx[i] bindings

Per AGENTS.md single-pathway rule: this script does NOT implement the upwind
stencil itself — it exercises the pipeline that any correct binding must match.

Run from the repo root:
  julia --project=. tests/conformance/discretization/rect_1d_advection_upwind_nonuniform_periodic/regenerate_golden.jl
=#

using JSON
using EarthSciDiscretizations: lower_stencil_to_replacement

# --------------------------------------------------------------------------
# Per-cell evaluator for non-uniform rules (handles index("dx", "$x"))
# --------------------------------------------------------------------------

function _resolve_idx(ie, cell_idx::Int)::Int
    ie isa AbstractString && String(ie) == "\$x" && return cell_idx
    ie isa AbstractDict || error("unsupported index expression: $ie")
    ia = ie["args"]
    String(ia[1]) == "\$x" || error("index LHS must be '\$x', got $(ia[1])")
    op2 = String(ie["op"])
    op2 == "+" && return cell_idx + Int(ia[2])
    op2 == "-" && return cell_idx - Int(ia[2])
    error("unsupported index op '$op2'")
end

function eval_nonuniform(node, u::Vector{Float64}, dx_arr::Vector{Float64},
                         cell_idx::Int, N::Int)::Float64
    node isa Number && return Float64(node)
    if node isa AbstractString
        s = String(node)
        s == "\$x" && return Float64(cell_idx)
        error("unresolved variable '$s'")
    end
    node isa AbstractDict || error("unexpected node type: $(typeof(node))")
    op   = String(node["op"])
    args = node["args"]
    if op == "index"
        arr_name = String(args[1])
        j = _resolve_idx(args[2], cell_idx)
        arr_name == "\$u" && return u[mod(j - 1, N) + 1]
        arr_name == "dx"  && return dx_arr[mod(j - 1, N) + 1]
        error("unknown array '$arr_name' in index op")
    end
    ev(a) = eval_nonuniform(a, u, dx_arr, cell_idx, N)
    op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
    op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
    op == "*" && return prod(ev(a) for a in args)
    op == "/" && return ev(args[1]) / ev(args[2])
    error("unsupported op '$op'")
end

function apply_rule_nonuniform(rule_path::String, rule_name::String,
                               u::Vector{Float64}, dx_arr::Vector{Float64})
    raw  = JSON.parsefile(rule_path)
    body = raw["discretizations"][rule_name]
    lw   = lower_stencil_to_replacement(body)
    repl = lw["replacement"]
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl
    N = length(u)
    return Float64[eval_nonuniform(expr, u, dx_arr, i, N) for i in 1:N]
end

# --------------------------------------------------------------------------
# Grid and IC helpers
# --------------------------------------------------------------------------

function stretched_dx(n::Int, x_start::Float64, x_end::Float64,
                      amplitude::Float64)::Vector{Float64}
    L = x_end - x_start
    h = L / n
    dx = Float64[h * (1.0 + amplitude * sin(2π * (i - 1) / n)) for i in 1:n]
    # Verify normalisation: sum(dx) must equal L to machine precision
    @assert abs(sum(dx) - L) < 1e-10 "stretched_dx: sum(dx)=$(sum(dx)) ≠ L=$L"
    return dx
end

function cell_centers(dx_arr::Vector{Float64}, x_start::Float64)::Vector{Float64}
    N = length(dx_arr)
    x_c = Vector{Float64}(undef, N)
    pos = x_start
    for i in 1:N
        x_c[i] = pos + dx_arr[i] / 2
        pos += dx_arr[i]
    end
    return x_c
end

function gaussian_field(x_c::Vector{Float64}, center::Float64,
                        sigma::Float64, scale::Float64)::Vector{Float64}
    amp = scale / (sigma * sqrt(2π))
    return Float64[amp * exp(-(xi - center)^2 / (2 * sigma^2)) for xi in x_c]
end

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

repo_root = abspath(joinpath(@__DIR__, "..", "..", "..", ".."))

fixtures_path = joinpath(@__DIR__, "fixtures.json")
FIXTURES = JSON.parsefile(fixtures_path)

rule_path = joinpath(repo_root, FIXTURES["rule_path"])

for fx in FIXTURES["fixtures"]
    g   = fx["grid"]
    ic  = fx["initial_condition"]
    N   = Int(g["n_cells"])
    x_start = Float64(g["x_start"])
    x_end   = Float64(g["x_end"])
    amp = Float64(g["stretching_amplitude"])

    dx_arr = stretched_dx(N, x_start, x_end, amp)
    x_c    = cell_centers(dx_arr, x_start)
    u      = gaussian_field(x_c, Float64(ic["center"]),
                            Float64(ic["sigma"]), Float64(ic["scale_factor"]))

    du_dx = apply_rule_nonuniform(rule_path, FIXTURES["rule"], u, dx_arr)

    golden = Dict(
        "_captured_by" => get(FIXTURES, "_captured_by", "esd-02z"),
        "field_u"      => u,
        "du_dx"        => du_dx,
        "dx"           => dx_arr,
    )

    out_path = joinpath(@__DIR__, "golden", "$(fx["name"]).json")
    open(out_path, "w") do io
        JSON.print(io, golden, 2)
    end
    println("Wrote golden: $out_path")
end
