@testsnippet DuoTopologyFAQSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON
    import EarthSciSerialization

    const _R = EarthSciSerialization.Relational

    # Serialization helpers — MUST match
    # tests/conformance/grids/duo/topology/regenerate_golden.jl byte-for-byte so
    # the golden comparison is apples-to-apples. Indexing is binding-neutral
    # 0-based; boundary sentinel -1 (Julia uses 1-based ids and 0 = no neighbor).
    _compact(x) = JSON.json(x)

    function _edges_ser(faq)
        pairs = [(faq.edges[1, e] - 1, faq.edges[2, e] - 1) for e in 1:size(faq.edges, 2)]
        return _R.canonical_index_set_json(pairs)
    end
    function _eof_ser(faq)
        Nc = size(faq.edges_on_face, 2)
        return _compact([[faq.edges_on_face[k, c] - 1 for k in 1:3] for c in 1:Nc])
    end
    function _cn_ser(faq)
        Nc = size(faq.cell_neighbors, 2)
        return _compact([
            [faq.cell_neighbors[k, c] == 0 ? -1 : faq.cell_neighbors[k, c] - 1 for k in 1:3]
            for c in 1:Nc
        ])
    end
    function _vf_ser(faq)
        return _compact([[f - 1 for f in faq.vertex_faces[v]] for v in 1:length(faq.vertex_faces)])
    end

    _build(level) = build_duo_grid(
        loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral", check = "strict"),
    )
end

@testitem "DUO primal topology value-invention FAQ — matches imperative + byte golden" setup = [DuoTopologyFAQSetup] tags = [:conformance, :grid, :duo, :topology, :faq] begin
    GOLDEN = JSON.parsefile(
        joinpath(@__DIR__, "..", "tests", "conformance", "grids", "duo", "topology", "golden.json"),
    )
    by_level = Dict(Int(l["level"]) => l for l in GOLDEN["levels"])

    for level in 0:2
        g = _build(level)
        faces = g.faces
        Nc = size(faces, 2)
        Nv = size(g.vertices, 2)
        faq = primal_topology_faq(faces, Nv)

        # --- counts ---
        @test size(faq.edges, 2) == n_edges(g)
        @test size(faq.cell_neighbors, 2) == n_cells(g)
        @test length(faq.vertex_faces) == n_vertices(g)

        # --- "match imperative output": cell_neighbors / vertex_faces are
        # order-invariant and must equal the imperative builders EXACTLY. ---
        @test faq.cell_neighbors == g.cell_neighbors
        @test faq.vertex_faces == g.vertex_faces

        # --- edges: same undirected set as the imperative builder; the FAQ
        # emits them in the canonical §5.5 sorted order (the cross-binding
        # byte-stable numbering that replaces Dict-insertion order). ---
        imp_edges = Set((g.edges[1, e], g.edges[2, e]) for e in 1:size(g.edges, 2))
        faq_edges = Set((faq.edges[1, e], faq.edges[2, e]) for e in 1:size(faq.edges, 2))
        @test imp_edges == faq_edges
        @test issorted([(faq.edges[1, e], faq.edges[2, e]) for e in 1:size(faq.edges, 2)])
        @test all(faq.edges[1, e] < faq.edges[2, e] for e in 1:size(faq.edges, 2))  # a < b

        # --- edges_on_face: each slot's dense id addresses the correct sorted
        # vertex pair, and the neighbour across that slot shares the same edge. ---
        for c in 1:Nc
            v1, v2, v3 = faces[1, c], faces[2, c], faces[3, c]
            for (k, (a, b)) in enumerate(((v2, v3), (v3, v1), (v1, v2)))
                e = faq.edges_on_face[k, c]
                want = a < b ? (a, b) : (b, a)
                @test (faq.edges[1, e], faq.edges[2, e]) == want
                c2 = faq.cell_neighbors[k, c]
                if c2 != 0
                    @test e in (faq.edges_on_face[1, c2], faq.edges_on_face[2, c2], faq.edges_on_face[3, c2])
                end
            end
        end

        # --- byte-identity: the canonical 0-based serialization is the
        # binding-neutral contract every binding's M3 engine must reproduce. ---
        gl = by_level[level]["serialized"]
        @test _edges_ser(faq) == gl["edges"]
        @test _eof_ser(faq) == gl["edges_on_face"]
        @test _cn_ser(faq) == gl["cell_neighbors"]
        @test _vf_ser(faq) == gl["vertex_faces"]
    end
end

@testitem "DUO primal topology FAQ anchors to the ESS §7.3 cross-binding golden" setup = [DuoTopologyFAQSetup] tags = [:conformance, :grid, :duo, :topology, :faq] begin
    # ESS pins the value-invention edge enumeration byte-identically across
    # Julia/Rust/Python in tests/conformance/determinism/manifest.json (fixture
    # `edge_enumeration`): the canonical 2-triangle mesh {1,2,3},{3,2,4} → the
    # unique undirected edge set, serialized "[[1,2],[1,3],[2,3],[2,4],[3,4]]"
    # (dense ids 0..4 canonical). Feeding the SAME mesh through the ESD bridge
    # must reproduce those exact bytes — proving ESD's `primal_topology_faq`
    # rides the same cross-binding determinism contract that guarantees the DUO
    # icosahedral output is byte-identical across the three bindings.
    faces = [1 3; 2 2; 3 4]        # (3, 2): column c is triangle c
    faq = primal_topology_faq(faces, 4)
    # edges keep the input's vertex labels (1-based here, as the ESS golden pins).
    pairs = [(faq.edges[1, e], faq.edges[2, e]) for e in 1:size(faq.edges, 2)]
    @test _R.canonical_index_set_json(pairs) == "[[1,2],[1,3],[2,3],[2,4],[3,4]]"
    # rank: Julia native 1-based; canonical 0-based = id − 1 ⇒ 0,1,2,3,4.
    dense0 = [faq.edges_on_face[k, c] - 1 for c in 1:2 for k in 1:3]
    @test sort(unique(dense0)) == [0, 1, 2, 3, 4]
    # The shared edge {2,3} collapses to one key incident on both triangles.
    @test faq.cell_neighbors[1, 1] == 2  # triangle 1 slot 1 = edge (2,3) → triangle 2
    @test faq.cell_neighbors[3, 2] == 1  # triangle 2 slot 3 = edge (2,3) → triangle 1
end

@testitem "DUO primal topology FAQ document is schema-valid" tags = [:conformance, :grid, :duo, :topology, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "duo", "rules", "primal_topology.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion.
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    # The two value-invention index sets are present and `derived` from the
    # aggregate nodes that mint them.
    model = file.models["DuoPrimalTopology"]
    @test haskey(model.index_sets, "edges")
    @test haskey(model.index_sets, "vertex_face_incidence")
end
