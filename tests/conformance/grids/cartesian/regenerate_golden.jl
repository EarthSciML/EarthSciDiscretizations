#!/usr/bin/env julia
# Regenerate golden accessor outputs for the cartesian conformance harness.
# Per docs/GRIDS_API.md §4.3, Julia is the reference binding for ULP ties.
#
# This script activates a temporary environment so it runs cleanly without
# touching the main Project.toml deps. JSON is pulled in just for IO.

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "..", "..", "..")); io = devnull)
    Pkg.add("JSON"; io = devnull)
end

using EarthSciDiscretizations
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN_DIR = joinpath(HERE, "golden")

# Build a CartesianGrid from a fixture opts dict. Handles the uniform
# (nx/ny/nz + extent) and non-uniform (edges) paths per GRIDS_API.md §2.
function build_grid(opts::AbstractDict)
    kw = Dict{Symbol, Any}()
    if haskey(opts, "edges")
        kw[:edges] = [Vector{Float64}(Float64.(e)) for e in opts["edges"]]
    else
        haskey(opts, "nx") && (kw[:nx] = Int(opts["nx"]))
        haskey(opts, "ny") && (kw[:ny] = Int(opts["ny"]))
        haskey(opts, "nz") && (kw[:nz] = Int(opts["nz"]))
        if haskey(opts, "extent")
            kw[:extent] = [(Float64(e[1]), Float64(e[2])) for e in opts["extent"]]
        end
    end
    kw[:ghosts] = Int(get(opts, "ghosts", 0))
    # dtype is fixed to Float64 for the conformance corpus (§2.2).
    return EarthSciDiscretizations.grids.cartesian(; kw...)
end

# Query points in the fixtures use 0-based indices (binding-neutral, matches
# Python/Rust/TS). Julia accessors are 1-based; convert at the boundary.
function qp_1idx(qp)
    return Tuple(Int(x) + 1 for x in qp)
end

# Emit neighbors as a sorted list-of-entries matching the Rust
# `CartesianNeighbor` shape (axis-major, low side before high side, 0-indexed).
# Julia's neighbors() returns a Dict keyed by (axis, side) with axis 1..N and
# Julia 1-based indices; we convert back to 0-indexed here.
function neighbors_0idx(grid, idx1)
    dict = neighbors(grid, idx1...)
    out = Vector{Dict{String, Any}}()
    N = ndims(grid)
    for d in 1:N
        for s in (-1, +1)
            haskey(dict, (d, s)) || continue
            nb1 = dict[(d, s)]
            push!(out, Dict{String, Any}(
                "axis" => d - 1,
                "side" => s,
                "index" => [Int(x) - 1 for x in nb1],
            ))
        end
    end
    return out
end

# Metric fields that apply at every `ndim`. See GRIDS_API.md §2 and
# src/grids/cartesian.jl for the full taxonomy.
const _COMMON_METRICS = ("volume", "jacobian")
const _AXIS_WIDTH = ("dx", "dy", "dz")
const _FACE_AREA = ("face_area_x", "face_area_y", "face_area_z")

function metric_scalar(grid, name::String, idx1)
    sym = Symbol(name)
    return Float64(metric_eval(grid, sym, idx1...))
end

function metric_tensor_g(grid, idx1)
    g = metric_eval(grid, :g, idx1...)
    return [[Float64(g[i][j]) for j in 1:length(g[i])] for i in 1:length(g)]
end

function build_output(fixture::AbstractDict)
    opts = fixture["opts"]
    grid = build_grid(opts)
    N = ndims(grid)

    qps = fixture["query_points"]
    n_qp = length(qps)

    centers = Vector{Vector{Float64}}(undef, n_qp)
    widths = Vector{Vector{Float64}}(undef, n_qp)
    volumes = Vector{Float64}(undef, n_qp)
    nbrs = Vector{Vector{Dict{String, Any}}}(undef, n_qp)
    gtensor = Vector{Vector{Vector{Float64}}}(undef, n_qp)

    # Only include the metric columns the fixture's ndim actually supports;
    # this keeps the golden shape self-describing and avoids encoding
    # axis-does-not-exist sentinels.
    metric_cols = Dict{String, Vector{Float64}}()
    for name in _COMMON_METRICS
        metric_cols[name] = Vector{Float64}(undef, n_qp)
    end
    for d in 1:N
        metric_cols[_AXIS_WIDTH[d]] = Vector{Float64}(undef, n_qp)
        metric_cols[_FACE_AREA[d]] = Vector{Float64}(undef, n_qp)
    end

    for (k, qp) in enumerate(qps)
        length(qp) == N ||
            error("fixture $(fixture["name"]): qp[$k] arity $(length(qp)) ≠ ndim $N")
        idx1 = qp_1idx(qp)
        c = cell_centers(grid, idx1...)
        w = cell_widths(grid, idx1...)
        centers[k] = [Float64(x) for x in c]
        widths[k] = [Float64(x) for x in w]
        volumes[k] = Float64(cell_volume(grid, idx1...))
        nbrs[k] = neighbors_0idx(grid, idx1)
        gtensor[k] = metric_tensor_g(grid, idx1)
        for name in keys(metric_cols)
            metric_cols[name][k] = metric_scalar(grid, name, idx1)
        end
    end

    out = Dict{String, Any}(
        "fixture" => fixture["name"],
        "family" => "cartesian",
        "ndim" => N,
        "opts" => opts,
        "n_cells" => prod([length(grid.widths[d]) for d in 1:N]),
        "query_points" => qps,
        "cell_centers" => centers,
        "cell_widths" => widths,
        "cell_volume" => volumes,
        "neighbors" => nbrs,
        "metric_g" => gtensor,
    )
    for (name, vec) in metric_cols
        out["metric_$name"] = vec
    end
    return out
end

function main()
    spec = JSON.parsefile(FIXTURES)
    isdir(GOLDEN_DIR) || mkpath(GOLDEN_DIR)
    for fx in spec["fixtures"]
        out = build_output(fx)
        path = joinpath(GOLDEN_DIR, "$(fx["name"]).json")
        open(path, "w") do io
            JSON.print(io, out, 2)
            write(io, "\n")
        end
        println("wrote $path")
    end
end

main()
