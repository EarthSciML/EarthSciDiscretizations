#!/usr/bin/env julia
# Regenerate golden accessor outputs for the arakawa conformance harness.
# Per docs/GRIDS_API.md §4.3, Julia is the reference binding for ULP ties.
#
# The script activates a temporary environment so it runs cleanly without
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

const LOCATIONS = ("cell_center", "u_edge", "v_edge", "corner")
const DIR_KEYS = ("W", "E", "S", "N")
const STAGGERS = ("A", "B", "C", "D", "E")
const METRICS = ("dx", "dy", "area")

# Julia's arakawa accessors are 1-indexed; the harness wire form is 0-indexed
# to match Python/Rust/TS. Convert at the boundary.

function loc_to_julia(loc::AbstractString)
    return loc == "cell_center" ? CellCenter :
        loc == "u_edge" ? UEdge :
        loc == "v_edge" ? VEdge :
        loc == "corner" ? Corner :
        throw(ArgumentError("unknown location $loc"))
end

function stagger_to_julia(s::AbstractString)
    return s == "A" ? ArakawaA :
        s == "B" ? ArakawaB :
        s == "C" ? ArakawaC :
        s == "D" ? ArakawaD :
        s == "E" ? ArakawaE :
        throw(ArgumentError("unknown stagger $s"))
end

function loc_name(loc::VarLocation)
    return loc === CellCenter ? "cell_center" :
        loc === UEdge ? "u_edge" :
        loc === VEdge ? "v_edge" :
        "corner"
end

# Accessor at any location, evaluated at 0-indexed (i, j). The Julia API is
# 1-based and splits the accessor across (cell_centers, u_face, v_face,
# corners); we pass separate per-location grids so u_face hits UEdge and
# v_face hits VEdge regardless of stagger.
function coord_0idx(grid_cc::ArakawaGrid, grid_c::ArakawaGrid,
        loc::VarLocation, i::Int, j::Int)
    if loc === CellCenter
        return cell_centers(grid_cc, i + 1, j + 1)
    elseif loc === Corner
        return corners(grid_cc, i + 1, j + 1)
    elseif loc === UEdge
        # Stagger C puts u on UEdge; u_face on g_c returns the u-edge coord.
        return u_face(grid_c, i + 1, j + 1)
    else  # VEdge
        return v_face(grid_c, i + 1, j + 1)
    end
end

# Neighbors at 0-indexed (loc, i, j); returns 0-indexed tuples or `nothing`.
function neighbors_0idx(grid::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
    w, e, s, n = neighbors(grid, loc, i + 1, j + 1)
    to_wire(t) = t === nothing ? nothing : [t[2] - 1, t[3] - 1]
    return (to_wire(w), to_wire(e), to_wire(s), to_wire(n))
end

function build_base(base_opts::AbstractDict)
    base_opts["family"] == "cartesian" || error("only cartesian base supported today")
    return CartesianBase(
        xlo = float(base_opts["xlo"]),
        xhi = float(base_opts["xhi"]),
        ylo = float(base_opts["ylo"]),
        yhi = float(base_opts["yhi"]),
        nx = Int(base_opts["nx"]),
        ny = Int(base_opts["ny"]),
    )
end

function build_output(fixture::AbstractDict)
    opts = fixture["opts"]
    base = build_base(opts["base"])
    dtype = opts["dtype"] == "float32" ? Float32 : Float64
    ghosts = Int(opts["ghosts"])
    staggers = fixture["staggers"]
    qps_by_loc = fixture["query_points"]

    # Geometry (coord/neighbors/metric) is stagger-independent — u_edge and
    # v_edge coordinates only depend on the base grid. We keep two reference
    # grids: one at stagger A for cell_center/corner queries (any stagger
    # works; A is simplest) and one at stagger C so u_face/v_face route to
    # UEdge/VEdge respectively.
    ref_a = EarthSciDiscretizations.grids.arakawa(
        base = base, stagger = :A, ghosts = ghosts, dtype = dtype,
    )
    ref_c = EarthSciDiscretizations.grids.arakawa(
        base = base, stagger = :C, ghosts = ghosts, dtype = dtype,
    )

    coords = Dict{String, Any}()
    nbrs = Dict{String, Any}()
    for lname in LOCATIONS
        loc = loc_to_julia(lname)
        qps = qps_by_loc[lname]
        xy = Vector{Vector{Float64}}(undef, length(qps))
        w_list = Vector{Any}(undef, length(qps))
        e_list = Vector{Any}(undef, length(qps))
        s_list = Vector{Any}(undef, length(qps))
        n_list = Vector{Any}(undef, length(qps))
        for (k, qp) in enumerate(qps)
            i, j = Int(qp[1]), Int(qp[2])
            x, y = coord_0idx(ref_a, ref_c, loc, i, j)
            xy[k] = [Float64(x), Float64(y)]
            w, e, s, n = neighbors_0idx(ref_a, loc, i, j)
            w_list[k] = w
            e_list[k] = e
            s_list[k] = s
            n_list[k] = n
        end
        coords[lname] = Dict("points" => qps, "xy" => xy)
        nbrs[lname] = Dict(
            "points" => qps,
            "W" => w_list, "E" => e_list, "S" => s_list, "N" => n_list,
        )
    end

    # Metric eval — cartesian metrics are stagger-independent and
    # constant-per-grid; sample at the cell_center query points.
    qps_cc = qps_by_loc["cell_center"]
    metric_values = Dict{String, Vector{Float64}}()
    for m in METRICS
        metric_values[m] = Vector{Float64}(undef, length(qps_cc))
    end
    for (k, qp) in enumerate(qps_cc)
        i, j = Int(qp[1]) + 1, Int(qp[2]) + 1
        metric_values["dx"][k] = Float64(metric_eval(ref_a, :dx, i, j))
        metric_values["dy"][k] = Float64(metric_eval(ref_a, :dy, i, j))
        metric_values["area"][k] = Float64(metric_eval(ref_a, :area, i, j))
    end

    # Per-stagger variable location + shape tables.
    per_stagger = Dict{String, Any}()
    for sname in STAGGERS
        sname in staggers || continue
        s = stagger_to_julia(sname)
        g_s = EarthSciDiscretizations.grids.arakawa(
            base = base, stagger = s, ghosts = ghosts, dtype = dtype,
        )
        h_loc, u_loc, v_loc = arakawa_variable_locations(s)
        per_stagger[sname] = Dict(
            "rotated" => (s === ArakawaE),
            "variable_locations" => Dict(
                "h" => loc_name(h_loc),
                "u" => loc_name(u_loc),
                "v" => loc_name(v_loc),
            ),
            "variable_shapes" => Dict(
                "h" => collect(variable_shape(g_s, :h)),
                "u" => collect(variable_shape(g_s, :u)),
                "v" => collect(variable_shape(g_s, :v)),
            ),
            "location_shapes" => Dict(
                "cell_center" => collect(arakawa_shape(g_s, CellCenter)),
                "u_edge" => collect(arakawa_shape(g_s, UEdge)),
                "v_edge" => collect(arakawa_shape(g_s, VEdge)),
                "corner" => collect(arakawa_shape(g_s, Corner)),
            ),
        )
    end

    return Dict(
        "fixture" => fixture["name"],
        "family" => "arakawa",
        "base_family" => "cartesian",
        "opts" => opts,
        "nx" => base.nx,
        "ny" => base.ny,
        "n_cells" => base.nx * base.ny,
        "coords" => coords,
        "neighbors" => nbrs,
        "metrics" => Dict(
            "points" => qps_cc,
            "dx" => metric_values["dx"],
            "dy" => metric_values["dy"],
            "area" => metric_values["area"],
        ),
        "staggers" => per_stagger,
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
