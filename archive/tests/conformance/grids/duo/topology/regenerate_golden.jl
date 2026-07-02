#!/usr/bin/env julia
# Regenerate the DUO primal-topology value-invention conformance golden.
#
# The golden pins the binding-neutral, byte-identical canonical serialization of
# the four primal-topology arrays (edges / edges_on_face / cell_neighbors /
# vertex_faces) that the landed M3 value-invention engine produces from each
# icosahedral level's `faces` connectivity, via `primal_topology_faq`
# (src/topology_faq.jl, declarative companion
# discretizations/grids/duo/rules/primal_topology.esm).
#
# Indexing is binding-neutral 0-based (Julia's 1-based ids and 0 "no neighbor"
# are converted at this boundary, exactly like tests/conformance/grids/duo/
# regenerate_golden.jl): cells/vertices/faces/edges are 0-based and the boundary
# sentinel is -1. Because the ESS relational engine is byte-identical across
# Julia/Rust/Python (CONFORMANCE_SPEC.md §5.5) and the `faces` connectivity is
# the same deterministic icosahedron in every binding, each binding's FAQ output
# serializes to these identical bytes.
#
# Usage:  julia --project=. tests/conformance/grids/duo/topology/regenerate_golden.jl

using EarthSciDiscretizations
using JSON
import EarthSciSerialization

const _R = EarthSciSerialization.Relational

"Compact JSON (no spaces) of a nested integer array — the canonical byte form."
_compact(x) = JSON.json(x)

"0-based undirected edge set as a canonical (sorted) index set, ESS-serialized."
function _edges_0based(faq)
    pairs = [(faq.edges[1, e] - 1, faq.edges[2, e] - 1) for e in 1:size(faq.edges, 2)]
    return pairs, _R.canonical_index_set_json(pairs)
end

"(3, Nc) dense-edge-id matrix → 0-based row-per-face nested array + compact bytes."
function _eof_0based(faq)
    Nc = size(faq.edges_on_face, 2)
    rows = [[faq.edges_on_face[k, c] - 1 for k in 1:3] for c in 1:Nc]
    return rows, _compact(rows)
end

"(3, Nc) neighbor matrix → 0-based row-per-face nested array (-1 boundary) + bytes."
function _cn_0based(faq)
    Nc = size(faq.cell_neighbors, 2)
    rows = [[faq.cell_neighbors[k, c] == 0 ? -1 : faq.cell_neighbors[k, c] - 1 for k in 1:3] for c in 1:Nc]
    return rows, _compact(rows)
end

"Ragged vertex→faces → 0-based list-per-vertex nested array + compact bytes."
function _vf_0based(faq)
    rows = [[f - 1 for f in faq.vertex_faces[v]] for v in 1:length(faq.vertex_faces)]
    return rows, _compact(rows)
end

function build_level(level::Int)
    g = build_duo_grid(
        loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral", check = "strict"),
    )
    faq = primal_topology_faq(g.faces, size(g.vertices, 2))
    _, edges_ser = _edges_0based(faq)
    _, eof_ser = _eof_0based(faq)
    _, cn_ser = _cn_0based(faq)
    _, vf_ser = _vf_0based(faq)
    return Dict(
        "level" => level,
        "n_cells" => n_cells(g),
        "n_vertices" => n_vertices(g),
        "n_edges" => n_edges(g),
        "serialized" => Dict(
            "edges" => edges_ser,
            "edges_on_face" => eof_ser,
            "cell_neighbors" => cn_ser,
            "vertex_faces" => vf_ser,
        ),
    )
end

function main()
    out = Dict(
        "description" =>
            "DUO primal-topology value-invention conformance golden (bead esd-heg.1 / D1a). " *
            "Pins the binding-neutral 0-based canonical serialization of edges / edges_on_face / " *
            "cell_neighbors / vertex_faces produced by the landed M3 relational engine " *
            "(primal_topology_faq -> EarthSciSerialization.Relational) for icosahedral levels 0/1/2. " *
            "edges is the ESS canonical undirected index set; the dense arrays are compact JSON in " *
            "the deterministic face/vertex order. Every binding's FAQ output MUST serialize to these " *
            "identical bytes (CONFORMANCE_SPEC.md §5.5).",
        "reference_binding" => "julia",
        "indexing" => "0-based binding-neutral; boundary sentinel -1; convert at the harness boundary",
        "rank_base_pin" => Dict("julia" => 1, "rust" => 0, "python" => 0),
        "levels" => [build_level(level) for level in 0:2],
    )
    dir = @__DIR__
    path = joinpath(dir, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    return println("wrote ", path)
end

main()
