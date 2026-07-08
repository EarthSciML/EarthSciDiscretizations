#!/usr/bin/env julia
# Regenerate the arakawa construction value-FAQ conformance golden.
#
# The golden pins the binding-neutral serialization of the arakawa staggered-grid
# construction the landed M1 elementwise FAQ produces from each conformance fixture,
# via `arakawa_construction_faq` (src/arakawa_faq.jl, declarative companion
# discretizations/grids/arakawa/rules/arakawa_construction.esm):
#   - the four staggered location coordinates (cell_center / u_edge / v_edge / corner),
#   - the cell-centred dx/dy/cell_volume,
#   - the row-major 0-sentinel neighbor maps and boundary masks, and
#   - the static A–E variable-location / shape table.
#
# Fixtures are SHARED with the accessor golden: this script reads the SAME
# tests/conformance/grids/arakawa/fixtures.json, so the FAQ construction and the
# accessor conformance are anchored to one fixture corpus.
#
# Sampling: the shared arakawa fixtures are large (the realistic base is 256x256), so
# — exactly like the accessor golden ../golden/*.json — the coordinate / neighbor /
# boundary arrays are pinned at the fixtures' `query_points` rather than materialized
# in full (the test still asserts the FULL arrays match the imperative builder to ULP
# in memory). The static A–E table is small and pinned in full per stagger.
#
# Indexing is binding-neutral: coordinate/metric floats are full-precision Float64
# (compared at the family relative tolerance); neighbor maps are 0-based with a -1
# sentinel; boundary masks are 0/1. Because EarthSciSerialization evaluates the
# elementwise arithmetic identically in every binding (CONFORMANCE_SPEC.md §5.5) and
# the grid parameters are the same, each binding's FAQ output reproduces these values.
#
# Usage:  julia --project=. tests/conformance/grids/arakawa/construction/regenerate_golden.jl

using EarthSciDiscretizations
const ESD = EarthSciDiscretizations
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "..", "fixtures.json")
const STAGGERS = ("A", "B", "C", "D", "E")

function build_base(base_opts::AbstractDict)
    base_opts["family"] == "cartesian" || error("only cartesian base supported today")
    return ESD.CartesianBase(
        xlo = Float64(base_opts["xlo"]), xhi = Float64(base_opts["xhi"]),
        ylo = Float64(base_opts["ylo"]), yhi = Float64(base_opts["yhi"]),
        nx = Int(base_opts["nx"]), ny = Int(base_opts["ny"]),
    )
end

loc_name(loc::VarLocation) =
    loc === CellCenter ? "cell_center" :
    loc === UEdge ? "u_edge" :
    loc === VEdge ? "v_edge" : "corner"

# Serialization helpers — MUST match test/test_arakawa_construction_faq.jl byte-for-byte
# so the golden comparison is apples-to-apples. Query-point sampled; coordinate floats
# are compact JSON strings, neighbor ids 0-based with -1 sentinel, masks 0/1.
_compact(x) = JSON.json(x)
function _coord_str(xs, ys, ni::Int, points)
    return _compact(
        [
            [
                    Float64(xs[Int(p[2]) * ni + Int(p[1]) + 1]),
                    Float64(ys[Int(p[2]) * ni + Int(p[1]) + 1]),
                ] for p in points
        ]
    )
end
function _nbr_str(arr, nx::Int, points)
    return _compact(
        [
            (id = arr[Int(p[2]) * nx + Int(p[1]) + 1]; id == 0 ? -1 : Int(id) - 1) for p in points
        ]
    )
end
function _mask_str(arr, nx::Int, points)
    return _compact([arr[Int(p[2]) * nx + Int(p[1]) + 1] ? 1 : 0 for p in points])
end

function build_fixture(f::AbstractDict)
    opts = f["opts"]
    base = build_base(opts["base"])
    dtype = get(opts, "dtype", "float64") == "float32" ? Float32 : Float64
    ghosts = Int(get(opts, "ghosts", 0))
    nx, ny = base.nx, base.ny
    qp = f["query_points"]

    # Geometry (coords / neighbors / boundary) is stagger-independent — build once at A.
    g0 = ESD.grids.arakawa(base = base, stagger = :A, ghosts = ghosts, dtype = dtype)
    faq0 = arakawa_construction_faq(g0)
    ccp = qp["cell_center"]

    # Per-stagger static table.
    per_stagger = Dict{String, Any}()
    for sname in STAGGERS
        sname in f["staggers"] || continue
        g = ESD.grids.arakawa(
            base = base, stagger = Symbol(sname), ghosts = ghosts, dtype = dtype,
        )
        faq = arakawa_construction_faq(g)
        per_stagger[sname] = Dict(
            "rotated" => faq.rotated,
            "variable_locations" => Dict(
                "h" => loc_name(faq.variable_locations.h),
                "u" => loc_name(faq.variable_locations.u),
                "v" => loc_name(faq.variable_locations.v),
            ),
            "variable_shapes" => Dict(
                "h" => collect(faq.variable_shapes.h),
                "u" => collect(faq.variable_shapes.u),
                "v" => collect(faq.variable_shapes.v),
            ),
            "location_shapes" => Dict(
                "cell_center" => collect(faq.location_shapes.cell_center),
                "u_edge" => collect(faq.location_shapes.u_edge),
                "v_edge" => collect(faq.location_shapes.v_edge),
                "corner" => collect(faq.location_shapes.corner),
            ),
        )
    end

    return Dict(
        "name" => f["name"],
        "nx" => nx,
        "ny" => ny,
        "n_cells" => nx * ny,
        "dx" => Float64(faq0.dx),
        "dy" => Float64(faq0.dy),
        "cell_volume" => Float64(faq0.cell_volume[1]),
        "coords" => Dict(
            "points" => qp,
            "cell_center" => _coord_str(faq0.coords.cell_center..., nx, qp["cell_center"]),
            "u_edge" => _coord_str(faq0.coords.u_edge..., nx + 1, qp["u_edge"]),
            "v_edge" => _coord_str(faq0.coords.v_edge..., nx, qp["v_edge"]),
            "corner" => _coord_str(faq0.coords.corner..., nx + 1, qp["corner"]),
        ),
        "neighbor_indices" => Dict(
            "points" => ccp,
            "x_minus" => _nbr_str(faq0.neighbor_minus[1], nx, ccp),
            "x_plus" => _nbr_str(faq0.neighbor_plus[1], nx, ccp),
            "y_minus" => _nbr_str(faq0.neighbor_minus[2], nx, ccp),
            "y_plus" => _nbr_str(faq0.neighbor_plus[2], nx, ccp),
        ),
        "boundary" => Dict(
            "points" => ccp,
            "x_lower" => _mask_str(faq0.boundary_lower[1], nx, ccp),
            "x_upper" => _mask_str(faq0.boundary_upper[1], nx, ccp),
            "y_lower" => _mask_str(faq0.boundary_lower[2], nx, ccp),
            "y_upper" => _mask_str(faq0.boundary_upper[2], nx, ccp),
        ),
        "staggers" => per_stagger,
    )
end

function main()
    spec = JSON.parsefile(FIXTURES)
    tol = get(spec, "tolerance", Dict("relative" => 1.0e-14))
    out = Dict(
        "description" =>
            "Arakawa construction value-FAQ conformance golden (bead esd-3we.4 / S4). " *
            "Pins the binding-neutral arakawa staggered-grid construction — the four " *
            "staggered location coordinates (cell_center/u_edge/v_edge/corner), the " *
            "cell-centred dx/dy/cell_volume, the row-major 0-sentinel neighbor maps and " *
            "boundary masks, and the static A–E variable-location/shape table — produced " *
            "by the landed M1 elementwise FAQ (arakawa_construction_faq -> " *
            "EarthSciSerialization AST evaluator) for the shared arakawa conformance " *
            "fixtures. Coordinate/metric floats are full-precision Float64 (sampled at the " *
            "fixtures' query_points, compared at the family relative tolerance); neighbor " *
            "maps are 0-based with a -1 sentinel; boundary masks are 0/1. Every binding's " *
            "FAQ output MUST reproduce these values (CONFORMANCE_SPEC.md §5.5).",
        "reference_binding" => "julia",
        "indexing" =>
            "coords/metrics full-precision Float64 (relative tolerance), sampled at " *
            "query_points; neighbor maps 0-based, sentinel -1; boundary masks 0/1; " *
            "location codes CellCenter/UEdge/VEdge/Corner; converted at the harness boundary",
        "tolerance" => tol,
        "fixtures" => [build_fixture(f) for f in spec["fixtures"]],
    )
    path = joinpath(HERE, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    return println("wrote ", path)
end

main()
