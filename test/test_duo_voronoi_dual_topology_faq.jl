@testsnippet DuoVoronoiDualTopologyFAQSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON

    # Serialization helpers — MUST match
    # tests/conformance/grids/duo/voronoi_dual_topology/regenerate_golden.jl
    # byte-for-byte so the golden comparison is apples-to-apples. Indexing is
    # binding-neutral 0-based; each cell's ring is RAGGED to n_edges_on_cell (the
    # angular order is the list order). cells_on_cell entries are 0-based primal
    # vertex ids; edges_on_cell entries are 0-based D1a canonical primal edge ids.
    _compact(x) = JSON.json(x)

    _neoc_ser(faq) = _compact(faq.n_edges_on_cell)
    function _coc_ser(faq)
        Nv = length(faq.n_edges_on_cell)
        return _compact([[faq.cells_on_cell[i, v] - 1 for i in 1:faq.n_edges_on_cell[v]] for v in 1:Nv])
    end
    function _eoc_ser(faq)
        Nv = length(faq.n_edges_on_cell)
        return _compact([[faq.edges_on_cell[i, v] - 1 for i in 1:faq.n_edges_on_cell[v]] for v in 1:Nv])
    end

    _build_primal(level) = build_duo_grid(
        loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral", check = "strict"),
    )
end

@testitem "DUO Voronoi-dual topology value-invention FAQ — matches imperative + byte golden" setup = [DuoVoronoiDualTopologyFAQSetup] tags = [:conformance, :grid, :duo, :mpas, :topology, :faq] begin
    GOLDEN = JSON.parsefile(
        joinpath(
            @__DIR__, "..", "tests", "conformance", "grids", "duo", "voronoi_dual_topology", "golden.json",
        ),
    )
    by_level = Dict(Int(l["level"]) => l for l in GOLDEN["levels"])

    # The MPAS Voronoi dual requires a subdivided primal (level ≥ 1).
    for level in 1:3
        g = _build_primal(level)
        Nv = size(g.vertices, 2)
        faq = voronoi_dual_topology_faq(g.vertices, g.faces, g.cell_cart; R = 6.371e6)

        # Imperative reference: the topology arrays the value-invention FAQ replaces.
        imp = EarthSciDiscretizations._duo_voronoi_dual(level; R = 6.371e6)

        # --- counts: one dual cell per primal vertex; closed icosahedral dual ---
        @test length(faq.n_edges_on_cell) == Nv
        @test size(faq.cells_on_cell) == (faq.max_edges, Nv)
        @test size(faq.edges_on_cell) == (faq.max_edges, Nv)
        @test faq.max_edges == 6
        @test count(==(5), faq.n_edges_on_cell) == 12      # the 12 pentagonal seeds
        @test count(==(6), faq.n_edges_on_cell) == Nv - 12 # the rest hexagonal

        # --- "match imperative _duo_voronoi_dual topology": the angular ordering
        # (cells_on_cell ring) and valence are byte-identical to the imperative
        # builder. cells_on_cell holds primal vertex ids, independent of the edge
        # numbering, so EXACT equality is the byte-identity contract. ---
        @test faq.n_edges_on_cell == imp.n_edges_on_cell
        @test faq.cells_on_cell == imp.cells_on_cell
        @test faq.max_edges == imp.max_edges

        # --- edges_on_cell: same angular slots as cells_on_cell, addressing the D1a
        # CANONICAL (sorted) primal edge {v, w}. The FAQ renumbers edges canonically
        # (replacing the imperative Dict-insertion numbering), so the contract is
        # CONSISTENCY — each edge id is the canonical id of the shared primal edge
        # {cell, neighbour} — not equality with the imperative ids. ---
        prim = primal_topology_faq(g.faces, Nv)
        Ne = size(prim.edges, 2)
        for v in 1:Nv
            kv = faq.n_edges_on_cell[v]
            for i in 1:kv
                w = faq.cells_on_cell[i, v]
                e = faq.edges_on_cell[i, v]
                @test 1 <= e <= Ne
                want = v < w ? (v, w) : (w, v)
                @test (prim.edges[1, e], prim.edges[2, e]) == want
            end
            # padding past the valence is exactly zero in both arrays
            for i in (kv + 1):faq.max_edges
                @test faq.cells_on_cell[i, v] == 0
                @test faq.edges_on_cell[i, v] == 0
            end
        end

        # --- each dual neighbour is symmetric: if w is a neighbour of v, v is a
        # neighbour of w (the dual cell graph is undirected). ---
        for v in 1:Nv
            for i in 1:faq.n_edges_on_cell[v]
                w = faq.cells_on_cell[i, v]
                @test v in (faq.cells_on_cell[j, w] for j in 1:faq.n_edges_on_cell[w])
            end
        end

        # --- byte-identity: the canonical 0-based serialization is the
        # binding-neutral contract every binding's M3 engine must reproduce. ---
        gl = by_level[level]["serialized"]
        @test _neoc_ser(faq) == gl["n_edges_on_cell"]
        @test _coc_ser(faq) == gl["cells_on_cell"]
        @test _eoc_ser(faq) == gl["edges_on_cell"]
    end
end

@testitem "DUO Voronoi-dual topology FAQ document is schema-valid" tags = [:conformance, :grid, :duo, :mpas, :topology, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(
        @__DIR__, "..", "discretizations", "grids", "duo", "rules", "voronoi_dual_topology.esm",
    )
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion; the explicit
    # `schema_errors` check mirrors the D1a primal-topology test (structural
    # "no defining equation" notices on aggregate-defined state vars are an
    # accepted limitation shared with the merged D1a document, not a schema error).
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    model = file.models["DuoVoronoiDualTopology"]
    # The angular-ordering KEY (the tangent-plane atan2 bearing) and its tangent
    # frame / incident-face inputs are present.
    @test haskey(model.index_sets, "cells")
    @test haskey(model.index_sets, "ring")
    @test haskey(model.variables, "angle")
    @test haskey(model.variables, "incident_face")
end
