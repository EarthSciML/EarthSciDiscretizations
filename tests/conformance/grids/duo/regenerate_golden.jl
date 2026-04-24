#!/usr/bin/env julia
# Regenerate golden accessor outputs for the duo conformance harness.
# Per docs/GRIDS_API.md §4.3, Julia is the reference binding for ULP ties.
#
# Activates a temporary env so it runs cleanly without touching the main
# Project.toml deps; JSON is pulled in just for IO.

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

const METRIC_NAMES = ("area", "lon", "lat", "x", "y", "z")
const BOUNDARY_SENTINEL = -1

# Julia's DuoGrid is 1-indexed internally and uses 0 to mean "no neighbor";
# the golden JSON is binding-neutral with 0-based cell indices and -1 for
# boundaries. Convert at the harness boundary so all four bindings compare
# byte-for-byte on the same numbers.
function neighbors_0idx(g, c1::Int)
    nb = EarthSciDiscretizations.neighbors(g, c1)
    return [n == 0 ? BOUNDARY_SENTINEL : n - 1 for n in nb]
end

function build_output(fixture::AbstractDict)
    opts = fixture["opts"]
    level = Int(opts["level"])
    R = float(opts["R"])
    ghosts = Int(opts["ghosts"])

    loader = DuoLoader(path = "builtin://icosahedral/$(level)")
    grid = build_duo_grid(; loader = loader, R = R, ghosts = ghosts)

    qps = Int.(fixture["query_points"])
    n = length(qps)

    centers = Vector{Dict{String, Float64}}(undef, n)
    nbrs = Vector{Vector{Int}}(undef, n)
    metrics = Dict(m => Vector{Float64}(undef, n) for m in METRIC_NAMES)

    for (k, c0) in enumerate(qps)
        c1 = c0 + 1
        lon, lat = EarthSciDiscretizations.cell_centers(grid, c1)
        centers[k] = Dict("lon" => lon, "lat" => lat)
        nbrs[k] = neighbors_0idx(grid, c1)
        for m in METRIC_NAMES
            metrics[m][k] = EarthSciDiscretizations.metric_eval(grid, Symbol(m), c1)
        end
    end

    return Dict(
        "fixture" => fixture["name"],
        "family" => "duo",
        "generator" => "builtin_icosahedral",
        "opts" => opts,
        "level" => level,
        "n_cells" => EarthSciDiscretizations.n_cells(grid),
        "n_vertices" => EarthSciDiscretizations.n_vertices(grid),
        "n_edges" => EarthSciDiscretizations.n_edges(grid),
        "query_points" => qps,
        "cell_centers" => centers,
        "neighbors" => nbrs,
        "metric_area" => metrics["area"],
        "metric_lon" => metrics["lon"],
        "metric_lat" => metrics["lat"],
        "metric_x" => metrics["x"],
        "metric_y" => metrics["y"],
        "metric_z" => metrics["z"],
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
