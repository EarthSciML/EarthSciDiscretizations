#!/usr/bin/env julia
# Regenerate the lat-lon construction value-FAQ conformance golden.
#
# The golden pins the binding-neutral serialization of the lat-lon grid construction
# arrays (per-cell lon/lat centers and widths, cell_volume, the 1-D latitude
# edges/centers, the closed-form curvilinear metric components, the metric Jacobian,
# the lone non-vanishing metric derivative, the periodic / nearest-center neighbor
# maps and the latitude boundary masks) that the landed M1 elementwise FAQ produces
# from each conformance fixture, via `latlon_construction_faq` (src/latlon_faq.jl,
# declarative companion discretizations/grids/latlon/rules/latlon_construction.esm).
#
# Fixtures are SHARED with the accessor golden: this script reads the SAME
# tests/conformance/grids/latlon/fixtures.json, so the FAQ construction and the
# accessor conformance are anchored to one fixture corpus (regular + reduced_gaussian).
#
# Indexing is binding-neutral: coordinate/metric floats are full-precision Float64
# (compared at the family's relative tolerance); neighbor maps are 0-based with a -1
# sentinel (Julia's 1-based flat ids and 0 = no-neighbor are converted at this
# boundary, exactly like tests/conformance/grids/cartesian/construction/
# regenerate_golden.jl); boundary masks are 0/1. Because EarthSciSerialization
# evaluates the elementwise arithmetic — including the sin/cos of the latitudes —
# identically in every binding (CONFORMANCE_SPEC.md §5.5) and the grid parameters are
# the same, each binding's FAQ output reproduces these values.
#
# Usage:  julia --project=. tests/conformance/grids/latlon/construction/regenerate_golden.jl

using EarthSciDiscretizations
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "..", "fixtures.json")

# Build a LatLonGrid from a fixture opts dict — the regular (nlon/nlat) and
# reduced_gaussian (nlon_per_row + supplied lat_edges) paths. MUST match
# test/test_latlon_construction_faq.jl byte-for-byte.
function build_grid(opts::AbstractDict)
    variant = Symbol(get(opts, "variant", "regular"))
    dtype = get(opts, "dtype", "float64") == "float32" ? Float32 : Float64
    kw = Dict{Symbol, Any}(:variant => variant, :dtype => dtype)
    kw[:R] = Float64(opts["R"])
    kw[:ghosts] = Int(get(opts, "ghosts", 0))
    haskey(opts, "pole_policy") && (kw[:pole_policy] = Symbol(opts["pole_policy"]))
    if variant === :regular
        kw[:nlon] = Int(opts["nlon"])
        kw[:nlat] = Int(opts["nlat"])
    else
        kw[:nlon_per_row] = Int[Int(x) for x in opts["nlon_per_row"]]
        haskey(opts, "lat_edges") &&
            (kw[:lat_edges] = Float64[Float64(x) for x in opts["lat_edges"]])
    end
    return EarthSciDiscretizations.grids.lat_lon(; kw...)
end

# Serialization helpers — MUST match test/test_latlon_construction_faq.jl byte-for-byte
# so the golden comparison is apples-to-apples. Each construction array is one compact
# JSON string of the flat ragged-row-major cell vector (or the 1-D latitude vector),
# the same dense byte form as tests/conformance/grids/cartesian/construction/golden.json.
"Compact JSON (no spaces) — the canonical byte form."
_compact(x) = JSON.json(x)
"Flat float array → compact full-precision Float64 string."
_floats(v) = _compact([Float64(x) for x in v])
"Tensor component [k, c...] over all cells → compact full-precision Float64 string."
_comp(a, idx...) = _compact([Float64(a[k, idx...]) for k in 1:size(a, 1)])
"Neighbor map → compact 0-based ids, -1 = no neighbor (0-sentinel converted)."
_nbr(v) = _compact([id == 0 ? -1 : Int(id) - 1 for id in v])
"Boundary mask → compact 0/1."
_mask(v) = _compact([b ? 1 : 0 for b in v])

function build_fixture(name::AbstractString, opts::AbstractDict)
    g = build_grid(opts)
    faq = latlon_construction_faq(g)
    return Dict(
        "name" => name,
        "opts" => opts,
        "variant" => String(g.variant),
        "nlat" => g.nlat,
        "nlon_per_row" => Int[x for x in g.nlon_per_row],
        "n_cells" => EarthSciDiscretizations.n_cells(g),
        "cell_centers_lon" => _floats(faq.cell_centers_lon),
        "cell_centers_lat" => _floats(faq.cell_centers_lat),
        "cell_widths_lon" => _floats(faq.cell_widths_lon),
        "cell_widths_lat" => _floats(faq.cell_widths_lat),
        "cell_volume" => _floats(faq.cell_volume),
        "lat_edges" => _floats(faq.lat_edges),
        "lat_centers" => _floats(faq.lat_centers),
        "metric_g_lonlon" => _comp(faq.metric_g, 1, 1),
        "metric_g_latlat" => _comp(faq.metric_g, 2, 2),
        "metric_ginv_lonlon" => _comp(faq.metric_ginv, 1, 1),
        "metric_jacobian" => _floats(faq.metric_jacobian),
        "dg_lonlon_dlat" => _comp(faq.metric_dgij_dxk, 1, 1, 2),
        "neighbor_lon_minus" => _nbr(faq.neighbor_lon_minus),
        "neighbor_lon_plus" => _nbr(faq.neighbor_lon_plus),
        "neighbor_lat_minus" => _nbr(faq.neighbor_lat_minus),
        "neighbor_lat_plus" => _nbr(faq.neighbor_lat_plus),
        "boundary_lat_lower" => _mask(faq.boundary_lat_lower),
        "boundary_lat_upper" => _mask(faq.boundary_lat_upper),
    )
end

function main()
    spec = JSON.parsefile(FIXTURES)
    tol = get(spec, "tolerance", Dict("relative" => 1.0e-14))
    out = Dict(
        "description" =>
            "Lat-lon construction value-FAQ conformance golden (bead esd-3we.3 / S2). " *
            "Pins the binding-neutral construction arrays — per-cell lon/lat centers and " *
            "widths, cell_volume, the 1-D latitude edges/centers, the closed-form " *
            "curvilinear metric (g_lonlon = R^2 cos^2 phi, g_latlat = R^2, ginv_lonlon, the " *
            "Jacobian R^2 |cos phi|, the lone derivative dg_lonlon/dphi = -2 R^2 cos phi sin " *
            "phi), the periodic / nearest-center neighbor maps and the latitude boundary " *
            "masks — produced by the landed M1 elementwise FAQ (latlon_construction_faq -> " *
            "EarthSciSerialization AST evaluator) for the shared lat-lon conformance " *
            "fixtures (regular + reduced_gaussian). Coordinate/metric floats are " *
            "full-precision Float64 compared at the family relative tolerance; neighbor maps " *
            "are 0-based with a -1 sentinel; boundary masks are 0/1. Every binding's FAQ " *
            "output MUST reproduce these values (CONFORMANCE_SPEC.md §5.5). Mirrors the " *
            "STRUCTURED-GRID FAQ TEMPLATE (cartesian esd-3we.1 / S1).",
        "reference_binding" => "julia",
        "reference_binding_notes" =>
            "Julia is the reference binding for the construction FAQ (it routes the same " *
            "arithmetic through EarthSciSerialization that every binding does). Distinct " *
            "from the accessor golden in ../golden/, whose reference_binding is python " *
            "pending the Julia lat_lon accessor landing (dsc-1ts).",
        "indexing" =>
            "coords/metrics full-precision Float64 (relative tolerance); neighbor maps " *
            "0-based, sentinel -1; boundary masks 0/1; flat ragged-row-major cell order; " *
            "converted at the harness boundary",
        "tolerance" => tol,
        "fixtures" => [build_fixture(f["name"], f["opts"]) for f in spec["fixtures"]],
    )
    path = joinpath(HERE, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    return println("wrote ", path)
end

main()
