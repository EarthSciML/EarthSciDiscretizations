@testsnippet ScvtTopologyLeafSetup begin
    using Test
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON

    # Serialization helpers — MUST match
    # tests/conformance/grids/mpas/scvt/topology_leaf/regenerate_golden.jl
    # byte-for-byte so the golden comparison is apples-to-apples. Indexing is
    # binding-neutral 0-based; each cell's ring is RAGGED to n_edges_on_cell.
    _compact(x) = JSON.json(x)
    _faces_ser(conn) = _compact([[conn.faces[k, t] - 1 for k in 1:3] for t in 1:conn.n_triangles])
    function _ring_ser(arr, neoc)
        Ncells = length(neoc)
        return _compact([[arr[i, c] - 1 for i in 1:neoc[c]] for c in 1:Ncells])
    end

    # The octahedron seed: 6 axis generators (unit magnitude; the leaf normalizes).
    # MUST match the regenerator's column order so the canonical golden applies.
    _octahedron() = Float64[
        1.0 -1.0  0.0  0.0  0.0  0.0
        0.0  0.0  1.0 -1.0  0.0  0.0
        0.0  0.0  0.0  0.0  1.0 -1.0
    ]

    _icosahedral_l1() = build_duo_grid(
        loader = (path = "builtin://icosahedral/1", reader = "builtin_icosahedral", check = "strict"),
    ).vertices

    # Canonical unordered face identity (a sorted triple of GENERATOR ids), the
    # order-independent fingerprint of the triangulation.
    _face_set(conn) = Set(Tuple(sort([conn.faces[1, t], conn.faces[2, t], conn.faces[3, t]])) for t in 1:conn.n_triangles)
end

@testitem "MPAS-SCVT topology leaf — octahedron: exact connectivity + byte golden" setup = [ScvtTopologyLeafSetup] tags = [:conformance, :grid, :mpas, :scvt, :topology, :leaf] begin
    GOLDEN = JSON.parsefile(
        joinpath(
            @__DIR__, "..", "tests", "conformance", "grids", "mpas", "scvt", "topology_leaf", "golden.json"
        )
    )
    seed = first(s for s in GOLDEN["seeds"] if s["name"] == "octahedron")
    gl = seed["serialized"]

    conn = scvt_voronoi_connectivity(_octahedron(); R = 1.0)

    # --- the smallest closed mesh: 6 generators → 2·6−4 = 8 triangles ---
    @test conn.n_triangles == 8
    @test size(conn.faces) == (3, 8)
    @test conn.n_edges_on_cell == [4, 4, 4, 4, 4, 4]   # every Voronoi cell a square
    @test conn.max_edges == 4

    # --- byte-identical canonical serialization (the determinism contract) ---
    @test _faces_ser(conn) == gl["faces"]
    @test _compact(conn.n_edges_on_cell) == gl["n_edges_on_cell"]
    @test _ring_ser(conn.cells_on_cell, conn.n_edges_on_cell) == gl["cells_on_cell"]
    @test _ring_ser(conn.edges_on_cell, conn.n_edges_on_cell) == gl["edges_on_cell"]
    @test _ring_ser(conn.vertices_on_cell, conn.n_edges_on_cell) == gl["vertices_on_cell"]

    # --- a cell neighbours every generator except its antipode (octahedron) ---
    # generators: 1=+x 2=−x 3=+y 4=−y 5=+z 6=−z; antipodal pairs (1,2)(3,4)(5,6).
    antipode = Dict(1 => 2, 2 => 1, 3 => 4, 4 => 3, 5 => 6, 6 => 5)
    for c in 1:6
        nbrs = Set(conn.cells_on_cell[1:conn.n_edges_on_cell[c], c])
        @test antipode[c] ∉ nbrs
        @test nbrs == Set(setdiff(1:6, [c, antipode[c]]))
    end
end

@testitem "MPAS-SCVT topology leaf — structural invariants + S2B vertex alignment" setup = [ScvtTopologyLeafSetup] tags = [:conformance, :grid, :mpas, :scvt, :topology, :leaf] begin
    for gens in (_octahedron(), _icosahedral_l1())
        conn = scvt_voronoi_connectivity(gens; R = 6.371e6)
        N = size(gens, 2)

        # Euler: a closed triangulated sphere of N generators has 2N−4 triangles.
        @test conn.n_triangles == 2 * N - 4

        # Every undirected Delaunay edge bounds exactly two triangles (closed mesh).
        edge_tri_count = Dict{Tuple{Int, Int}, Int}()
        for t in 1:conn.n_triangles
            a, b, c = conn.faces[1, t], conn.faces[2, t], conn.faces[3, t]
            for (u, v) in ((a, b), (b, c), (a, c))
                key = u < v ? (u, v) : (v, u)
                edge_tri_count[key] = get(edge_tri_count, key, 0) + 1
            end
        end
        @test all(==(2), values(edge_tri_count))

        # Canonical per-triangle form: smallest index first, CCW-from-outside
        # (outward normal · vertex > 0), and the list lexicographically sorted.
        tri_list = [(conn.faces[1, t], conn.faces[2, t], conn.faces[3, t]) for t in 1:conn.n_triangles]
        @test issorted(tri_list)
        for (a, b, c) in tri_list
            @test a < b && a < c            # smallest first
            pa = gens[:, a]; pb = gens[:, b]; pc = gens[:, c]
            n = ESD._cross3(Tuple(pb .- pa), Tuple(pc .- pa))
            @test ESD._dot3(n, Tuple(pa)) > 0   # CCW as seen from outside
        end

        # Dual-graph adjacency is symmetric (undirected Voronoi/Delaunay edges).
        for c in 1:N, i in 1:conn.n_edges_on_cell[c]
            w = conn.cells_on_cell[i, c]
            @test c in (conn.cells_on_cell[j, w] for j in 1:conn.n_edges_on_cell[w])
        end

        # vertices_on_cell ↔ cells_on_cell S2B alignment: slot i is the triangle
        # between neighbours i and (i%k)+1.
        for c in 1:N
            k = conn.n_edges_on_cell[c]
            for i in 1:k
                t = conn.vertices_on_cell[i, c]
                want = Set((c, conn.cells_on_cell[i, c], conn.cells_on_cell[(i % k) + 1, c]))
                @test Set(conn.faces[:, t]) == want
            end
        end

        # Tolerance geometry: every Voronoi vertex (circumcentre) lies on the
        # sphere of radius R (the geometry tolerance contract, NOT byte-pinned).
        for t in 1:conn.n_triangles
            nrm = sqrt(conn.circumcenters[1, t]^2 + conn.circumcenters[2, t]^2 + conn.circumcenters[3, t]^2)
            @test isapprox(nrm, 6.371e6; rtol = 1.0e-12, atol = 1.0e-3)
        end
    end
end

@testitem "MPAS-SCVT topology leaf — determinism: permutation-invariant connectivity" setup = [ScvtTopologyLeafSetup] tags = [:conformance, :grid, :mpas, :scvt, :topology, :leaf, :determinism] begin
    # The integer connectivity is a pure function of the generator SET, not the
    # input order: relabel the generators by a permutation and the triangulation
    # — as a set of generator-identity triples — is unchanged. This is the
    # cross-binding determinism property the contract requires (the canonical
    # ordering removes any hull-construction-order dependence).
    for (gens, perm) in (
            (_octahedron(), [4, 1, 6, 2, 5, 3]),
            (_icosahedral_l1(), vcat(22:42, 1:21)),   # a fixed rotation of the 42 ids
        )
        base = scvt_voronoi_connectivity(gens; R = 6.371e6)

        permuted_gens = gens[:, perm]
        pconn = scvt_voronoi_connectivity(permuted_gens; R = 6.371e6)

        # Map each permuted face's new ids back to original ids, as a sorted
        # triple, and compare the SET to the original triangulation's face set.
        mapped = Set(
            Tuple(sort([perm[pconn.faces[1, t]], perm[pconn.faces[2, t]], perm[pconn.faces[3, t]]]))
                for t in 1:pconn.n_triangles
        )
        @test mapped == _face_set(base)

        # Re-running on identical input is byte-identical (trivial determinism).
        @test scvt_voronoi_connectivity(gens; R = 6.371e6).faces == base.faces
    end
end

@testitem "MPAS-SCVT topology leaf — icosahedral rho≡1 regression vs imperative _duo_voronoi_dual" setup = [ScvtTopologyLeafSetup] tags = [:conformance, :grid, :mpas, :scvt, :topology, :leaf, :regression] begin
    GOLDEN = JSON.parsefile(
        joinpath(
            @__DIR__, "..", "tests", "conformance", "grids", "mpas", "scvt", "topology_leaf", "golden.json"
        )
    )
    seed = first(s for s in GOLDEN["seeds"] if s["name"] == "icosahedral_level1")
    gl = seed["serialized"]

    g = build_duo_grid(
        loader = (path = "builtin://icosahedral/1", reader = "builtin_icosahedral", check = "strict")
    )
    R = 6.371e6
    conn = scvt_voronoi_connectivity(g.vertices; R = R)
    imp = ESD._duo_voronoi_dual(1; R = R)
    Nv = size(g.vertices, 2)

    # Uniform-density SCVT (= CVT) on the icosahedral seed IS the quasi-uniform
    # icosahedral-dual MPAS mesh: 42 cells, 80 dual vertices, 12 pentagons.
    @test Nv == 42
    @test conn.n_triangles == 80
    @test count(==(5), conn.n_edges_on_cell) == 12
    @test count(==(6), conn.n_edges_on_cell) == 30

    # The headline rho≡1 regression: the leaf's `cells_on_cell` (generator-id
    # ring) is BYTE-IDENTICAL to the landed imperative `_duo_voronoi_dual` — the
    # convex-hull spherical Delaunay recovers exactly the DUO-dual topology.
    @test conn.n_edges_on_cell == imp.n_edges_on_cell
    @test conn.cells_on_cell == imp.cells_on_cell

    # And byte-identical to the cross-binding canonical golden.
    @test _faces_ser(conn) == gl["faces"]
    @test _ring_ser(conn.cells_on_cell, conn.n_edges_on_cell) == gl["cells_on_cell"]
    @test _ring_ser(conn.edges_on_cell, conn.n_edges_on_cell) == gl["edges_on_cell"]
    @test _ring_ser(conn.vertices_on_cell, conn.n_edges_on_cell) == gl["vertices_on_cell"]
end

@testitem "MPAS-SCVT topology leaf — degenerate / malformed inputs rejected" setup = [ScvtTopologyLeafSetup] tags = [:grid, :mpas, :scvt, :topology, :leaf] begin
    # Fewer than 4 generators cannot form a closed mesh.
    @test_throws ArgumentError scvt_spherical_delaunay(Float64[1 0 0; 0 1 0; 0 0 1]; R = 1.0)
    # All-coplanar generators (a great circle) are degenerate / non-enclosing.
    coplanar = Float64[1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0; 0.0 0.0 0.0 0.0]
    @test_throws Exception scvt_spherical_delaunay(coplanar; R = 1.0)
    # Bad radius / shape.
    @test_throws DomainError scvt_spherical_delaunay(_octahedron(); R = -1.0)
    @test_throws ArgumentError scvt_spherical_delaunay(Float64[1.0 0.0; 0.0 1.0]; R = 1.0)
end

@testitem "MPAS-SCVT topology leaf declaration document is schema-valid" tags = [:conformance, :grid, :mpas, :scvt, :topology, :leaf, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "mpas", "scvt", "topology_leaf.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion (mirrors the D1
    # background-quadrature + voronoi_dual_topology tests; structural notices on
    # parameter-only / aggregate-defined state vars are an accepted limitation).
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    model = file.models["MpasScvtTopologyLeaf"]
    @test haskey(model.index_sets, "triangles")
    @test haskey(model.variables, "triangle_cell")
    @test haskey(model.variables, "circumcenter")
end
