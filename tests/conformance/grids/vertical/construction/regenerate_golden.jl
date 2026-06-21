#!/usr/bin/env julia
# Regenerate the vertical construction value-FAQ conformance golden.
#
# The golden pins the binding-neutral serialization of the vertical column
# construction arrays (interface levels, cell centers/widths, cell_volume, the
# named layer metrics dz/z/sigma/pressure/ak/bk, the 0-sentinel plus/minus-k
# neighbor maps, the boundary masks, and the pressure_coefficients) that the landed
# M1 elementwise FAQ produces from each conformance fixture, via
# `vertical_construction_faq` (src/vertical_faq.jl, declarative companion
# discretizations/grids/vertical/rules/vertical_construction.esm).
#
# Fixtures are SHARED with the existing `.esm` conformance: this script reads the
# SAME tests/conformance/grids/vertical/fixtures.json, so the FAQ construction and
# the to_esm() conformance are anchored to one fixture corpus.
#
# Indexing is binding-neutral: level/center/width/metric floats are full-precision
# Float64 (compared at the family's relative tolerance — 0.0, i.e. strict byte
# equality, since the vertical family is pure rationals + a single multiply-add);
# neighbor maps are 0-based with a -1 sentinel (Julia's 1-based ids and 0 =
# no-neighbor are converted at this boundary, exactly like the cartesian
# construction golden); boundary masks are 0/1. Because EarthSciSerialization
# evaluates the elementwise arithmetic identically in every binding
# (CONFORMANCE_SPEC.md §5.5) and the grid parameters are the same, each binding's FAQ
# output reproduces these values.
#
# Usage:  julia --project=. tests/conformance/grids/vertical/construction/regenerate_golden.jl

using EarthSciDiscretizations
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "..", "fixtures.json")

# Build a VerticalGrid from a fixture opts dict — the uniform (coordinate + nz),
# supplied-levels (coordinate + levels) and hybrid (coordinate + ak/bk/p0) paths.
function build_grid(opts::AbstractDict)
    kw = Dict{Symbol, Any}(:coordinate => Symbol(opts["coordinate"]))
    haskey(opts, "nz") && (kw[:nz] = Int(opts["nz"]))
    haskey(opts, "levels") && (kw[:levels] = Float64.(opts["levels"]))
    haskey(opts, "ak") && (kw[:ak] = Float64.(opts["ak"]))
    haskey(opts, "bk") && (kw[:bk] = Float64.(opts["bk"]))
    haskey(opts, "p0") && (kw[:p0] = Float64(opts["p0"]))
    haskey(opts, "ghosts") && (kw[:ghosts] = Int(opts["ghosts"]))
    return EarthSciDiscretizations.grids.vertical(; kw...)
end

# Serialization helpers — MUST match test/test_vertical_construction_faq.jl
# byte-for-byte so the golden comparison is apples-to-apples. Each construction
# array is one compact JSON string (no spaces), the same dense byte form as the
# cartesian construction golden.
"Compact JSON (no spaces) — the canonical byte form."
_compact(x) = JSON.json(x)
"Float vector → compact full-precision Float64 string."
_floats(v) = _compact(Float64.(v))
"Neighbor map → compact 0-based linear ids, -1 = no neighbor."
_nbr(v) = _compact([id == 0 ? -1 : Int(id) - 1 for id in v])
"Boundary mask → compact 0/1."
_mask(v) = _compact([b ? 1 : 0 for b in v])

function build_fixture(name::AbstractString, opts::AbstractDict)
    g = build_grid(opts)
    faq = vertical_construction_faq(g)
    rec = Dict{String, Any}(
        "name" => name,
        "opts" => opts,
        "coordinate" => String(g.coordinate),
        "n_cells" => n_cells(g),
        "n_vertices" => n_vertices(g),
        "n_edges" => n_edges(g),
        "levels" => _floats(faq.levels),
        "centers" => _floats(faq.centers),
        "widths" => _floats(faq.widths),
        "cell_volume" => _floats(faq.cell_volume),
        "metric_dz" => _floats(faq.metric_dz),
        "metric_z" => _floats(faq.metric_z),
        "neighbor_minus" => _nbr(faq.neighbor_minus),
        "neighbor_plus" => _nbr(faq.neighbor_plus),
        "boundary_lower" => _mask(faq.boundary_lower),
        "boundary_upper" => _mask(faq.boundary_upper),
        "p0" => Float64(faq.p0),
    )
    faq.metric_sigma === nothing || (rec["metric_sigma"] = _floats(faq.metric_sigma))
    faq.metric_pressure === nothing || (rec["metric_pressure"] = _floats(faq.metric_pressure))
    faq.metric_ak === nothing || (rec["metric_ak"] = _floats(faq.metric_ak))
    faq.metric_bk === nothing || (rec["metric_bk"] = _floats(faq.metric_bk))
    isempty(faq.ak) || (rec["ak"] = _floats(faq.ak))
    isempty(faq.bk) || (rec["bk"] = _floats(faq.bk))
    return rec
end

function main()
    spec = JSON.parsefile(FIXTURES)
    tol = get(spec, "tolerance", Dict("relative" => 0.0))
    out = Dict(
        "description" =>
            "Vertical construction value-FAQ conformance golden (bead esd-3we.2 / S2). " *
            "Pins the binding-neutral construction arrays — interface levels, cell " *
            "centers/widths, cell_volume, the named layer metrics (dz/z/sigma/pressure/ak/bk), " *
            "the 0-sentinel plus/minus-k neighbor maps, the boundary masks and the " *
            "pressure_coefficients (ak/bk/p0) — produced by the landed M1 elementwise FAQ " *
            "(vertical_construction_faq -> EarthSciSerialization AST evaluator) for the shared " *
            "vertical conformance fixtures. Level/center/width/metric floats are full-precision " *
            "Float64 compared at the family relative tolerance (0.0 = strict byte equality); " *
            "neighbor maps are 0-based with a -1 sentinel; boundary masks are 0/1. Every " *
            "binding's FAQ output MUST reproduce these values (CONFORMANCE_SPEC.md §5.5). Reuses " *
            "the STRUCTURED-GRID FAQ TEMPLATE established by cartesian (esd-3we.1 / S1).",
        "reference_binding" => "julia",
        "indexing" =>
            "levels/centers/widths/metrics full-precision Float64 (relative tolerance); neighbor " *
            "maps 0-based, sentinel -1; boundary masks 0/1; converted at the harness boundary",
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
