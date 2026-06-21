#!/usr/bin/env julia
# Regenerate the MPAS-SCVT spherical-topology LEAF conformance golden
# (epic esd-e5m, bead esd-e5m.3 / D3).
#
# The golden pins the binding-neutral, byte-identical canonical serialization of
# the DETERMINISTIC INTEGER CONNECTIVITY the spherical-topology leaf emits from a
# set of generators — the spherical Delaunay `faces` and the dual Voronoi
# connectivity (`n_edges_on_cell` / `cells_on_cell` / `edges_on_cell` /
# `vertices_on_cell`) produced by `scvt_voronoi_connectivity`
# (src/grids/mpas_scvt_topology.jl): the irreducible convex-hull spherical
# Delaunay (canonical executor the s2bindings.rs S2B FFI, bead s2b-s7b) composed
# with the landed declarative dual-topology FAQ (`voronoi_dual_topology_faq`,
# esd-heg.2). Declarative companion + contract:
# discretizations/grids/mpas/scvt/topology_leaf.esm and
# discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md.
#
# Two canonical seeds, both NON-degenerate simplicial polytopes (no four
# cospherically-coplanar generators), so the Float64 convex hull is exact and the
# emitted integer connectivity is the unique, order-independent topology every
# binding (including the exact-predicate S2B FFI) MUST reproduce byte-for-byte:
#
#   * octahedron  — 6 axis generators → 8 Delaunay triangles (2·6−4), every
#                   Voronoi cell a 4-neighbour square. The smallest closed mesh.
#   * icosahedral level 1 — the 42 DUO primal vertices → 80 triangles
#                   (2·42−4), 12 pentagonal + 30 hexagonal cells. This seed makes
#                   the rho≡1 SCVT regression concrete: the leaf's `cells_on_cell`
#                   is byte-identical to the imperative `_duo_voronoi_dual`
#                   (asserted in test/test_mpas_scvt_topology_leaf.jl).
#
# Indexing is binding-neutral 0-based (Julia's 1-based ids converted at this
# boundary, exactly like the voronoi_dual_topology golden): `faces`,
# `cells_on_cell` (neighbour generator ids), `vertices_on_cell` (dual Voronoi
# vertex = triangle ids) are 0-based; `edges_on_cell` is the canonical (sorted)
# primal-edge numbering, 0-based. The per-cell rings are serialized RAGGED to
# n_edges_on_cell (angular order; the padding-to-max_edges convention is
# reconstructible from n_edges_on_cell and not pinned).
#
# The circumcentre GEOMETRY is NOT pinned here — it rides the TOLERANCE contract
# (its unit-norm-times-R invariant is asserted to tolerance in the test), not the
# byte-identity determinism contract that governs the integer connectivity.
#
# Usage:  julia --project=. tests/conformance/grids/mpas/scvt/topology_leaf/regenerate_golden.jl

using EarthSciDiscretizations
using JSON

const R_EARTH = 6.371e6

"Compact JSON (no spaces) of a nested integer array — the canonical byte form."
_compact(x) = JSON.json(x)

"""
The octahedron seed: 6 generators at the ±x/±y/±z axes (unit magnitude; the leaf
normalizes). Column order is fixed so the canonical connectivity is reproducible.
"""
function octahedron_generators()
    return Float64[
        1.0 -1.0  0.0  0.0  0.0  0.0
        0.0  0.0  1.0 -1.0  0.0  0.0
        0.0  0.0  0.0  0.0  1.0 -1.0
    ]
end

"0-based list of the canonical Delaunay triangles (each a 3-tuple of cell ids)."
function _faces_ser(conn)
    nt = conn.n_triangles
    return _compact([[conn.faces[k, t] - 1 for k in 1:3] for t in 1:nt])
end

"Ragged 0-based per-cell ring serializer (cells_on_cell / vertices_on_cell)."
function _ring_ser(arr, n_edges_on_cell)
    Ncells = length(n_edges_on_cell)
    return _compact([[arr[i, c] - 1 for i in 1:n_edges_on_cell[c]] for c in 1:Ncells])
end

function build_seed(name::String, generators::AbstractMatrix; R::Real = R_EARTH)
    conn = scvt_voronoi_connectivity(generators; R = R)
    return Dict(
        "name" => name,
        "n_cells" => size(generators, 2),
        "n_triangles" => conn.n_triangles,
        "max_edges" => conn.max_edges,
        "serialized" => Dict(
            "faces" => _faces_ser(conn),
            "n_edges_on_cell" => _compact(conn.n_edges_on_cell),
            "cells_on_cell" => _ring_ser(conn.cells_on_cell, conn.n_edges_on_cell),
            "edges_on_cell" => _ring_ser(conn.edges_on_cell, conn.n_edges_on_cell),
            "vertices_on_cell" => _ring_ser(conn.vertices_on_cell, conn.n_edges_on_cell),
        ),
    )
end

function icosahedral_l1_generators()
    g = build_duo_grid(
        loader = (path = "builtin://icosahedral/1", reader = "builtin_icosahedral", check = "strict"),
    )
    return g.vertices
end

function main()
    out = Dict(
        "description" =>
            "MPAS-SCVT spherical-topology LEAF conformance golden (epic esd-e5m, bead esd-e5m.3 / " *
            "D3). Pins the binding-neutral 0-based canonical serialization of the DETERMINISTIC " *
            "integer connectivity the leaf emits — the spherical Delaunay `faces` and the dual " *
            "Voronoi `n_edges_on_cell` / `cells_on_cell` / `edges_on_cell` / `vertices_on_cell` — " *
            "produced by `scvt_voronoi_connectivity` (the convex-hull spherical Delaunay, canonical " *
            "executor the s2bindings.rs S2B FFI s2b-s7b, composed with the declarative " *
            "voronoi_dual_topology_faq esd-heg.2). Two NON-degenerate seeds: the octahedron (6→8) " *
            "and icosahedral level 1 (42→80). Every binding's leaf output MUST serialize to these " *
            "identical bytes (determinism contract, CONFORMANCE_SPEC.md §5.5 / §5.7; " *
            "discretizations/grids/mpas/scvt/TOPOLOGY_LEAF_CONTRACT.md). The circumcentre geometry " *
            "rides the tolerance contract and is asserted to tolerance in " *
            "test/test_mpas_scvt_topology_leaf.jl, not pinned here. The icosahedral-L1 " *
            "`cells_on_cell` is byte-identical to the imperative `_duo_voronoi_dual` (the rho≡1 " *
            "SCVT regression), also asserted there.",
        "reference_binding" => "julia",
        "indexing" =>
            "0-based binding-neutral; faces / cells_on_cell = generator (cell) ids; " *
            "vertices_on_cell = dual Voronoi vertex (Delaunay triangle) ids; edges_on_cell = " *
            "canonical (sorted) primal-edge ids; rings ragged to n_edges_on_cell; convert at the " *
            "harness boundary",
        "rank_base_pin" => Dict("julia" => 1, "rust" => 0, "python" => 0),
        "seeds" => [
            build_seed("octahedron", octahedron_generators()),
            build_seed("icosahedral_level1", icosahedral_l1_generators()),
        ],
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
