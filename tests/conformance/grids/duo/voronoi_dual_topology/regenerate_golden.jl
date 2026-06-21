#!/usr/bin/env julia
# Regenerate the DUO→MPAS Voronoi-dual TOPOLOGY value-invention conformance golden.
#
# The golden pins the binding-neutral, byte-identical canonical serialization of
# the three dual-topology arrays — n_edges_on_cell / cells_on_cell / edges_on_cell
# — that the value-invention FAQ produces from each icosahedral level's DUO primal
# mesh, via `voronoi_dual_topology_faq` (src/topology_faq.jl, declarative companion
# discretizations/grids/duo/rules/voronoi_dual_topology.esm). This is the topology
# half of the imperative `_duo_voronoi_dual` (src/grids/mpas.jl, steps 3-4),
# expressed declaratively (bead esd-heg.2 / D1b).
#
# Indexing is binding-neutral 0-based (Julia's 1-based ids are converted at this
# boundary, exactly like tests/conformance/grids/duo/topology/regenerate_golden.jl):
# dual cells (= primal vertices) and the neighbour vertex ids in cells_on_cell are
# 0-based; the edge ids in edges_on_cell are the D1a CANONICAL (sorted) primal edge
# numbering, 0-based. Each cell's ring is serialized RAGGED — the active entries in
# angular order, length n_edges_on_cell — so the padding-to-max_edges convention is
# not pinned (it is reconstructible from n_edges_on_cell). The MPAS dual requires a
# subdivided primal, so the levels are 1/2/3 (n_cells = 42 / 162 / 642).
#
# Because the ESS relational engine (the bridge-vertex intersection / edge-id rank)
# is byte-identical across Julia/Rust/Python (CONFORMANCE_SPEC.md §5.5) and the
# angular `atan2` sort rides the geometry pipeline's fixed evaluation order, each
# binding's FAQ output serializes to these identical bytes.
#
# Usage:  julia --project=. tests/conformance/grids/duo/voronoi_dual_topology/regenerate_golden.jl

using EarthSciDiscretizations
using JSON

"Compact JSON (no spaces) of a nested integer array — the canonical byte form."
_compact(x) = JSON.json(x)

"Valence per dual cell (= primal-vertex degree); identity, no rebasing."
function _neoc_ser(faq)
    return _compact(faq.n_edges_on_cell)
end

"Ragged 0-based neighbour-vertex ring per cell, in angular order + compact bytes."
function _coc_ser(faq)
    Nv = length(faq.n_edges_on_cell)
    rows = [[faq.cells_on_cell[i, v] - 1 for i in 1:faq.n_edges_on_cell[v]] for v in 1:Nv]
    return _compact(rows)
end

"Ragged 0-based canonical-edge-id ring per cell, in angular order + compact bytes."
function _eoc_ser(faq)
    Nv = length(faq.n_edges_on_cell)
    rows = [[faq.edges_on_cell[i, v] - 1 for i in 1:faq.n_edges_on_cell[v]] for v in 1:Nv]
    return _compact(rows)
end

function build_level(level::Int)
    g = build_duo_grid(
        loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral", check = "strict"),
    )
    faq = voronoi_dual_topology_faq(g.vertices, g.faces, g.cell_cart; R = 6.371e6)
    return Dict(
        "level" => level,
        "n_cells" => size(g.vertices, 2),
        "n_dual_vertices" => size(g.faces, 2),
        "max_edges" => faq.max_edges,
        "serialized" => Dict(
            "n_edges_on_cell" => _neoc_ser(faq),
            "cells_on_cell" => _coc_ser(faq),
            "edges_on_cell" => _eoc_ser(faq),
        ),
    )
end

function main()
    out = Dict(
        "description" =>
            "DUO→MPAS Voronoi-dual TOPOLOGY value-invention conformance golden (bead esd-heg.2 / " *
            "D1b). Pins the binding-neutral 0-based canonical serialization of n_edges_on_cell / " *
            "cells_on_cell / edges_on_cell produced by the value-invention FAQ " *
            "(voronoi_dual_topology_faq -> EarthSciSerialization.Relational + the tangent-plane " *
            "atan2 angular sort) for icosahedral dual levels 1/2/3. cells_on_cell holds the " *
            "angularly-ordered bridge-vertex ring (0-based primal vertex ids); edges_on_cell holds " *
            "the D1a canonical (sorted) primal edge ids of those rings (0-based); both are ragged " *
            "to n_edges_on_cell. Every binding's FAQ output MUST serialize to these identical bytes " *
            "(CONFORMANCE_SPEC.md §5.5). The imperative `_duo_voronoi_dual` parity " *
            "(n_edges_on_cell / cells_on_cell byte-identical) is asserted in " *
            "test/test_duo_voronoi_dual_topology_faq.jl.",
        "reference_binding" => "julia",
        "indexing" =>
            "0-based binding-neutral; cells_on_cell = primal vertex ids; edges_on_cell = D1a " *
            "canonical (sorted) primal edge ids; rings ragged to n_edges_on_cell; convert at the " *
            "harness boundary",
        "rank_base_pin" => Dict("julia" => 1, "rust" => 0, "python" => 0),
        "levels" => [build_level(level) for level in 1:3],
    )
    dir = @__DIR__
    path = joinpath(dir, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote ", path)
end

main()
