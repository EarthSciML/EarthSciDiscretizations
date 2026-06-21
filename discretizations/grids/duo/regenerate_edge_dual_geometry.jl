#!/usr/bin/env julia
# Re-pin the D2b edge+dual GEOMETRY golden
# (faq/fixtures/canonical/edge_dual_geometry_level1.json) to the ESS CANONICAL
# (sorted) primal-edge order — the cross-binding-stable numbering that
# `primal_topology_faq` (D1a value-invention) now produces in `build_duo_grid`,
# aligned with the already-canonical D1a topology golden and the Rust/Python
# front-door parity landed in ess-3lj.2 (F2).
#
# This is a PURE PERMUTATION of the existing golden values: only the five
# edge-indexed arrays (dc_edge, dv_edge, cell_distance, lon_edge, lat_edge) are
# reordered from the legacy imperative `_build_connectivity` Dict-insertion order
# to the canonical order. Per-face (circ_x/y/z) and per-cell (lon/lat/area_cell)
# arrays are untouched. The Float64 VALUES are byte-identical — verified four ways
# below — so the cross-binding byte/value contract is preserved; only the edge
# index assignment changes (the explicit goal of the value-invention migration,
# bead esd-ohd / W1).
#
# The output is emitted as `JSON.print(io, nt, 2)` over a NamedTuple with
# alphabetically-sorted keys and NO trailing newline — byte-identical to the
# original generator's format, so the diff is exactly the five reordered arrays.
#
# Run from the repo root:
#     julia --project=. discretizations/grids/duo/regenerate_edge_dual_geometry.jl

using EarthSciDiscretizations
const ESD = EarthSciDiscretizations
using JSON

const GOLDEN = joinpath(
    @__DIR__, "faq", "fixtures", "canonical", "edge_dual_geometry_level1.json",
)
const EDGE_KEYS = ["dc_edge", "dv_edge", "cell_distance", "lon_edge", "lat_edge"]

bits(x) = reinterpret(Int64, Float64(x))

function main()
    orig_text = read(GOLDEN, String)
    old = JSON.parsefile(GOLDEN)
    level = Int(old["level"])

    # Canonical grid (build_duo_grid now wires primal_topology_faq → canonical edges).
    g = build_duo_grid(loader = (path = "builtin://icosahedral/$(level)", reader = "builtin_icosahedral"))
    ec = g.edges                                  # canonical (2, Ne)
    Ne = size(ec, 2)

    # Legacy imperative edge order — the order the existing golden was pinned to.
    _, F = ESD.duo_subdivide_faq(Float64, level)  # same faces build_duo_grid uses
    ei, _ = ESD._build_connectivity(F)            # imperative (2, Ne)
    @assert size(ei, 2) == Ne "imperative/canonical edge counts differ"

    # perm[ce] = imperative index of the physical edge that canonical index ce names.
    imp_idx = Dict((ei[1, i], ei[2, i]) => i for i in 1:Ne)
    perm = [imp_idx[(ec[1, ce], ec[2, ce])] for ce in 1:Ne]
    @assert sort(perm) == collect(1:Ne) "perm is not a bijection — edge sets differ"

    new = Dict{String, Any}()
    for (k, v) in old
        new[k] = (k in EDGE_KEYS) ? [v[perm[ce]] for ce in 1:Ne] : v
    end

    # --- Verification 1: each reordered array is a pure permutation (multiset kept).
    for k in EDGE_KEYS
        @assert sort(Float64.(new[k])) == sort(Float64.(old[k])) "value multiset changed for $k"
    end
    # --- Verification 2: reordered dc_edge == edge_length(g) (canonical), bit-for-bit.
    el = edge_length(g)
    @assert all(bits(new["dc_edge"][e]) == bits(el[e]) for e in 1:Ne) "dc_edge reorder != edge_length(g)"
    # --- Verification 3: reordered cell_distance == cell_distance(g), bit-for-bit.
    cd = cell_distance(g)
    @assert all(bits(new["cell_distance"][e]) == bits(cd[e]) for e in 1:Ne) "cell_distance reorder != cell_distance(g)"

    # Emit with alphabetically-sorted keys, indent 2, no trailing newline —
    # byte-identical format to the original generator.
    ks = sort(collect(keys(new)))
    nt = NamedTuple{Tuple(Symbol.(ks))}(Tuple(new[k] for k in ks))
    io = IOBuffer()
    JSON.print(io, nt, 2)
    out = String(take!(io))

    # --- Verification 4: only the five edge arrays' element ORDER changed; every
    # other byte (keys, formatting, per-face/per-cell values) is identical.
    new_unchanged = JSON.parse(out)
    for k in keys(old)
        k in EDGE_KEYS && continue
        if old[k] isa AbstractVector
            @assert all(bits(new_unchanged[k][i]) == bits(old[k][i]) for i in eachindex(old[k])) "untouched array $k changed"
        end
    end
    @assert length(out) == length(orig_text) "length changed beyond a pure reorder ($(length(orig_text)) → $(length(out)))"

    write(GOLDEN, out)
    println("re-pinned $(GOLDEN)")
    println("  5 edge arrays canonically reordered (pure permutation, multiset preserved);")
    println("  dc_edge == edge_length(g) and cell_distance == cell_distance(g) bit-for-bit;")
    println("  per-face/per-cell arrays + formatting byte-identical; file length unchanged ($(length(out)) bytes).")
    return
end

main()
