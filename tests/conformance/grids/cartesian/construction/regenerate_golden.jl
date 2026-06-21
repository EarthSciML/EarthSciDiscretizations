#!/usr/bin/env julia
# Regenerate the cartesian construction value-FAQ conformance golden.
#
# The golden pins the binding-neutral serialization of the cartesian grid
# construction arrays (per-axis edges/centers/widths, cell_volume, the identity
# metric Jacobian, the 0-sentinel neighbor maps and the boundary masks) that the
# landed M1 elementwise FAQ produces from each conformance fixture, via
# `cartesian_construction_faq` (src/cartesian_faq.jl, declarative companion
# discretizations/grids/cartesian/rules/cartesian_construction.esm).
#
# Fixtures are SHARED with the accessor golden: this script reads the SAME
# tests/conformance/grids/cartesian/fixtures.json, so the FAQ construction and the
# accessor conformance are anchored to one fixture corpus.
#
# Indexing is binding-neutral: coordinate/metric floats are full-precision Float64
# (compared at the family's relative tolerance); neighbor maps are 0-based with a
# -1 sentinel (Julia's 1-based linear ids and 0 = no-neighbor are converted at this
# boundary, exactly like tests/conformance/grids/cartesian/regenerate_golden.jl);
# boundary masks are 0/1. Because EarthSciSerialization evaluates the elementwise
# arithmetic identically in every binding (CONFORMANCE_SPEC.md §5.5) and the grid
# parameters are the same, each binding's FAQ output reproduces these values.
#
# Usage:  julia --project=. tests/conformance/grids/cartesian/construction/regenerate_golden.jl

using EarthSciDiscretizations
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "..", "fixtures.json")

# Build a CartesianGrid from a fixture opts dict — the uniform (nx/ny/nz + extent)
# and non-uniform (edges) paths, identical to the accessor regenerate script.
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
    return EarthSciDiscretizations.grids.cartesian(; kw...)
end

# Serialization helpers — MUST match test/test_cartesian_construction_faq.jl
# byte-for-byte so the golden comparison is apples-to-apples. Each construction
# array family is one compact JSON string (no spaces) of the [axis][...] nesting,
# the same dense byte form as tests/conformance/grids/duo/topology/golden.json.
"Compact JSON (no spaces) — the canonical byte form."
_compact(x) = JSON.json(x)
"Per-axis float arrays → compact [[axis 1 …],[axis 2 …],…] (full-precision Float64)."
_axis_floats(faq, field, N) = _compact([Float64.(getfield(faq, field)[d]) for d in 1:N])
"Per-axis neighbor maps → compact [[axis …]] of 0-based linear ids, -1 = no neighbor."
_axis_nbr(faq, field, N) =
    _compact([[id == 0 ? -1 : Int(id) - 1 for id in getfield(faq, field)[d]] for d in 1:N])
"Per-axis boundary masks → compact [[axis …]] of 0/1."
_axis_mask(faq, field, N) =
    _compact([[b ? 1 : 0 for b in getfield(faq, field)[d]] for d in 1:N])

function build_fixture(name::AbstractString, opts::AbstractDict)
    g = build_grid(opts)
    faq = cartesian_construction_faq(g)
    N = ndims(g)
    return Dict(
        "name" => name,
        "opts" => opts,
        "ndim" => N,
        "n" => collect(g.n),
        "n_cells" => prod(g.n),
        "edges" => _axis_floats(faq, :edges, N),
        "centers" => _axis_floats(faq, :centers, N),
        "widths" => _axis_floats(faq, :widths, N),
        "cell_volume" => _compact(Float64.(faq.cell_volume)),
        "metric_jacobian" => _compact(Float64.(faq.metric_jacobian)),
        "neighbor_minus" => _axis_nbr(faq, :neighbor_minus, N),
        "neighbor_plus" => _axis_nbr(faq, :neighbor_plus, N),
        "boundary_lower" => _axis_mask(faq, :boundary_lower, N),
        "boundary_upper" => _axis_mask(faq, :boundary_upper, N),
    )
end

function main()
    spec = JSON.parsefile(FIXTURES)
    tol = get(spec, "tolerance", Dict("relative" => 1e-14))
    out = Dict(
        "description" =>
            "Cartesian construction value-FAQ conformance golden (bead esd-3we.1 / S1). " *
            "Pins the binding-neutral construction arrays — per-axis edges/centers/widths, " *
            "cell_volume, the identity-metric Jacobian, the 0-sentinel neighbor maps and the " *
            "boundary masks — produced by the landed M1 elementwise FAQ " *
            "(cartesian_construction_faq -> EarthSciSerialization AST evaluator) for the shared " *
            "cartesian conformance fixtures. Coordinate/metric floats are full-precision Float64 " *
            "compared at the family relative tolerance; neighbor maps are 0-based with a -1 " *
            "sentinel; boundary masks are 0/1. Every binding's FAQ output MUST reproduce these " *
            "values (CONFORMANCE_SPEC.md §5.5). The STRUCTURED-GRID FAQ TEMPLATE that S2-S4 reuse.",
        "reference_binding" => "julia",
        "indexing" =>
            "coords/metrics full-precision Float64 (relative tolerance); neighbor maps 0-based, " *
            "sentinel -1; boundary masks 0/1; converted at the harness boundary",
        "tolerance" => tol,
        "fixtures" => [build_fixture(f["name"], f["opts"]) for f in spec["fixtures"]],
    )
    path = joinpath(HERE, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote ", path)
end

main()
