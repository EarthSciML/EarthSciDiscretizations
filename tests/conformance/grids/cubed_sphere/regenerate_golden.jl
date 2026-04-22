#!/usr/bin/env julia
# Regenerate golden accessor outputs for the cubed_sphere conformance harness.
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

const DIRECTIONS = (West, East, South, North)
const DIR_KEYS = ("W", "E", "S", "N")

# 0-indexed neighbor (p', i', j') from 0-indexed (p, i, j) on an Nc x Nc panel.
# Mirrors the Python/Rust/TS transform_ghost_index logic; cross-panel neighbors
# are comparable byte-for-byte across bindings.
function neighbor_0idx(p::Int, i::Int, j::Int, dir, Nc::Int)
    if dir == West && i > 0
        return (p, i - 1, j)
    elseif dir == East && i < Nc - 1
        return (p, i + 1, j)
    elseif dir == South && j > 0
        return (p, i, j - 1)
    elseif dir == North && j < Nc - 1
        return (p, i, j + 1)
    end
    nb = PANEL_CONNECTIVITY[p + 1][dir]
    nb_panel = nb.neighbor_panel - 1
    nb_edge = nb.neighbor_edge
    along = (dir == West || dir == East) ? j : i
    if nb.reverse_index
        along = Nc - 1 - along
    end
    if nb_edge == West
        return (nb_panel, 0, along)
    elseif nb_edge == East
        return (nb_panel, Nc - 1, along)
    elseif nb_edge == South
        return (nb_panel, along, 0)
    else
        return (nb_panel, along, Nc - 1)
    end
end

function cell_center_0idx(grid, p::Int, i::Int, j::Int)
    ξc = grid.ξ_centers[i + 1]
    ηc = grid.η_centers[j + 1]
    return gnomonic_to_lonlat(ξc, ηc, p + 1)
end

function metric_values_0idx(grid, p::Int, i::Int, j::Int)
    ξc = grid.ξ_centers[i + 1]
    ηc = grid.η_centers[j + 1]
    J, gxx, gyy, gxy = gnomonic_metric(ξc, ηc, grid.R)
    det = gxx * gyy - gxy * gxy
    return (
        J = J,
        g_xixi = gxx,
        g_etaeta = gyy,
        g_xieta = gxy,
        ginv_xixi = gyy / det,
        ginv_etaeta = gxx / det,
        ginv_xieta = -gxy / det,
    )
end

function cell_area_0idx(grid, p::Int, i::Int, j::Int)
    ξw = grid.ξ_edges[i + 1]
    ξe = grid.ξ_edges[i + 2]
    ηs = grid.η_edges[j + 1]
    ηn = grid.η_edges[j + 2]
    return compute_cell_area((ξw, ξe), (ηs, ηn), grid.R, p + 1)
end

function build_output(fixture::AbstractDict)
    opts = fixture["opts"]
    Nc = Int(opts["Nc"])
    R = float(opts["R"])
    ghosts = Int(opts["ghosts"])
    grid = CubedSphereGrid(Nc; R = R, Ng = ghosts)

    qps = fixture["query_points"]
    n = length(qps)
    centers = Vector{Dict{String,Float64}}(undef, n)
    nbrs = Dict(k => Vector{Vector{Int}}(undef, n) for k in DIR_KEYS)
    metric_names = ("J", "g_xixi", "g_etaeta", "g_xieta",
                    "ginv_xixi", "ginv_etaeta", "ginv_xieta")
    metrics = Dict(m => Vector{Float64}(undef, n) for m in metric_names)
    areas = Vector{Float64}(undef, n)

    for (k, qp) in enumerate(qps)
        p, i, j = Int(qp[1]), Int(qp[2]), Int(qp[3])
        lon, lat = cell_center_0idx(grid, p, i, j)
        centers[k] = Dict("lon" => lon, "lat" => lat)
        for (dkey, dir) in zip(DIR_KEYS, DIRECTIONS)
            np, ni, nj = neighbor_0idx(p, i, j, dir, Nc)
            nbrs[dkey][k] = [np, ni, nj]
        end
        m = metric_values_0idx(grid, p, i, j)
        metrics["J"][k] = m.J
        metrics["g_xixi"][k] = m.g_xixi
        metrics["g_etaeta"][k] = m.g_etaeta
        metrics["g_xieta"][k] = m.g_xieta
        metrics["ginv_xixi"][k] = m.ginv_xixi
        metrics["ginv_etaeta"][k] = m.ginv_etaeta
        metrics["ginv_xieta"][k] = m.ginv_xieta
        areas[k] = cell_area_0idx(grid, p, i, j)
    end

    return Dict(
        "fixture" => fixture["name"],
        "family" => "cubed_sphere",
        "generator" => "gnomonic_c6",
        "opts" => opts,
        "n_cells" => 6 * Nc * Nc,
        "query_points" => qps,
        "cell_centers" => centers,
        "neighbors_W" => nbrs["W"],
        "neighbors_E" => nbrs["E"],
        "neighbors_S" => nbrs["S"],
        "neighbors_N" => nbrs["N"],
        "metric_J" => metrics["J"],
        "metric_g_xixi" => metrics["g_xixi"],
        "metric_g_etaeta" => metrics["g_etaeta"],
        "metric_g_xieta" => metrics["g_xieta"],
        "metric_ginv_xixi" => metrics["ginv_xixi"],
        "metric_ginv_etaeta" => metrics["ginv_etaeta"],
        "metric_ginv_xieta" => metrics["ginv_xieta"],
        "area" => areas,
    )
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
