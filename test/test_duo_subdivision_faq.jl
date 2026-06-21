# DUO icosahedral subdivision expressed as a value-invention FAQ pass (esd-heg.3).
#
# One level of refinement is the RFC `semiring-faq-unified-ir` §7.3 value-invention
# pattern applied to a triangular mesh: per face, mint a midpoint vertex keyed on
# the canonical Skolem of its (sorted) endpoint pair, `distinct`-dedup the shared
# midpoints, `rank` them to dense ids, give each the normalized-average coordinate,
# join every face to its three edge-midpoints and emit four sub-faces. The declarative
# spec is `discretizations/grids/duo/faq/subdivide_refine.esm`; this test drives the
# *landed ESS engine primitives* (no shadow evaluator — AGENTS.md single-pathway) to
# materialize the pass and proves it reproduces the imperative `_subdivide_icosahedron`
# step exactly: same vertex coordinates, same triangles position-for-position (winding
# preserved). The only difference from the imperative arrays is the deterministic,
# cross-binding-canonical labeling of the invented vertices (§5.7 sorted order — the
# imperative `push!`/`Dict` order is not cross-binding stable, which is why the epic
# moves to FAQ).

@testitem "duo subdivision via value-invention FAQ matches imperative" tags = [:grid, :duo, :faq, :subdivision] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const REL = EarthSciSerialization.Relational

    # Midpoint coordinate AST: (v_a[d] + v_b[d]) / sqrt(Σ_d (v_a[d]+v_b[d])²).
    # Squares are products (`*`), NOT `^`, so the float ops match the imperative
    # `_midpoint!` bit-for-bit (Julia `mx^2` lowers to `mx*mx`; a float `^` would not).
    _sum(d) = Dict{String,Any}("op" => "+", "args" => ["a$d", "b$d"])
    _sq(d) = Dict{String,Any}("op" => "*", "args" => [_sum(d), _sum(d)])
    _norm() = Dict{String,Any}("op" => "sqrt",
        "args" => [Dict{String,Any}("op" => "+", "args" => [_sq(1), _sq(2), _sq(3)])])
    _component(d) = Dict{String,Any}("op" => "/", "args" => [_sum(d), _norm()])

    # One FAQ refine pass over a mesh (V 3×Nv, F 3×Nc) using the ESS engine:
    #   skolem_edge + distinct + rank (value invention) and eval_coeff (geometry FAQ).
    function faq_refine_once(V::Matrix{Float64}, F::Matrix{Int})
        Nv = size(V, 2)
        Nc = size(F, 2)
        # Directed local edges (a,b),(b,c),(c,a) of every face → canonical Skolem keys.
        pairs = Tuple{Int,Int}[]
        for f in 1:Nc
            a, b, c = F[1, f], F[2, f], F[3, f]
            push!(pairs, REL.skolem_edge(a, b))
            push!(pairs, REL.skolem_edge(b, c))
            push!(pairs, REL.skolem_edge(c, a))
        end
        # distinct + dense rank over the unique undirected edges (= midpoint vertices).
        rk = REL.rank(pairs; base = 1)          # rk.order = canonical sorted edges; rk.id = dense id
        edges = rk.order
        Ne = length(edges)

        # New vertices: base verts then the invented midpoints, in canonical order.
        Vnew = Matrix{Float64}(undef, 3, Nv + Ne)
        Vnew[:, 1:Nv] = V
        for (i, (lo, hi)) in enumerate(edges)
            b = Dict{String,Float64}(
                "a1" => V[1, lo], "a2" => V[2, lo], "a3" => V[3, lo],
                "b1" => V[1, hi], "b2" => V[2, hi], "b3" => V[3, hi],
            )
            Vnew[1, Nv + i] = ESD.eval_coeff(_component(1), b)
            Vnew[2, Nv + i] = ESD.eval_coeff(_component(2), b)
            Vnew[3, Nv + i] = ESD.eval_coeff(_component(3), b)
        end

        midvid(x, y) = Nv + rk.id[REL.skolem_edge(x, y)]   # base=1 ⇒ id ∈ 1..Ne

        # Per-face join to its 3 midpoints; emit 4 sub-faces in the imperative order.
        Fnew = Matrix{Int}(undef, 3, 4 * Nc)
        for f in 1:Nc
            a, b, c = F[1, f], F[2, f], F[3, f]
            ab = midvid(a, b); bc = midvid(b, c); ca = midvid(c, a)
            tris = ((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca))
            for (t, tri) in enumerate(tris)
                col = 4 * (f - 1) + t
                Fnew[1, col] = tri[1]; Fnew[2, col] = tri[2]; Fnew[3, col] = tri[3]
            end
        end
        return Vnew, Fnew, edges
    end

    coordset(V) = Set(Tuple(V[:, j]) for j in 1:size(V, 2))
    facetris(V, F) = [(Tuple(V[:, F[1, k]]), Tuple(V[:, F[2, k]]), Tuple(V[:, F[3, k]])) for k in 1:size(F, 2)]

    # A single FAQ pass applied to level n must reproduce the imperative level n+1
    # mesh exactly (n = 0 and n = 1).
    for n in 0:1
        Vn, Fn = ESD._subdivide_icosahedron(Float64, n)
        Vexp, Fexp = ESD._subdivide_icosahedron(Float64, n + 1)
        Vfaq, Ffaq, edges = faq_refine_once(Vn, Fn)

        # Counts: Nv' = Nv + Ne, Nc' = 4·Nc; Ne = 3·Nc/2 (closed triangular mesh).
        @test size(Vfaq, 2) == size(Vexp, 2)
        @test size(Ffaq, 2) == size(Fexp, 2)
        @test length(edges) == 3 * size(Fn, 2) ÷ 2

        # (A) Same vertex coordinate set, bit-for-bit (midpoint geometry FAQ ≡ _midpoint!).
        @test coordset(Vfaq) == coordset(Vexp)

        # (B) Same triangles, position-for-position with winding preserved — the FAQ
        #     emits sub-faces in the imperative order, so only midpoint *labels* differ.
        @test facetris(Vfaq, Ffaq) == facetris(Vexp, Fexp)

        # (C) Cross-binding determinism: the invented edge set is in §5.7 sorted order
        #     (never first-seen), independent of face iteration / winding.
        @test issorted(edges)
        @test allunique(edges)
        @test all(e -> e[1] < e[2], edges)   # canonical undirected (min,max)
    end

    # (D) Orientation/permutation independence of the value invention: shuffling and
    #     reversing the input faces must not change the invented edge set (the
    #     determinism contract that makes the labeling cross-binding stable).
    V0, F0 = ESD._subdivide_icosahedron(Float64, 0)
    _, _, edges0 = faq_refine_once(V0, F0)
    Fperm = F0[:, reverse(1:size(F0, 2))]            # permuted faces
    Frev = vcat(F0[3:3, :], F0[2:2, :], F0[1:1, :])  # reversed winding per face
    _, _, edges_perm = faq_refine_once(V0, Fperm)
    _, _, edges_rev = faq_refine_once(V0, Frev)
    @test edges0 == edges_perm
    @test edges0 == edges_rev
end

@testitem "duo subdivision FAQ: declarative doc + cross-binding canonical byte contract" tags = [:grid, :duo, :faq, :subdivision, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const REL = EarthSciSerialization.Relational
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    FAQ_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo", "faq")

    # --- The declarative value-invention FAQ document exists and expresses the chain.
    doc = JSON.parsefile(joinpath(FAQ_DIR, "subdivide_refine.esm"); dicttype = Dict{String, Any})
    model = doc["models"]["IcosahedralRefine"]
    isets = model["index_sets"]
    @test isets["edges"]["kind"] == "derived"
    @test isets["edges"]["from_faq"] == "edge_set"
    eqs = model["equations"]
    # edge_set: the value-invention aggregate — distinct + Skolem key, addressable by id.
    edge_set = only(filter(e -> get(e["rhs"], "id", nothing) == "edge_set", eqs))["rhs"]
    @test edge_set["distinct"] == true
    @test edge_set["key"]["op"] == "skolem"
    # rank node and the four-child emission must both be present.
    @test any(e -> e["rhs"]["expr"] isa Dict && get(e["rhs"]["expr"], "op", "") == "rank", eqs)
    @test any(e -> e["lhs"] isa Dict && e["lhs"]["args"][1] == "child_vert", eqs)

    # --- Cross-binding canonical byte contract for the invented edge/midpoint set.
    golden = JSON.parsefile(
        joinpath(FAQ_DIR, "fixtures", "canonical", "edge_set_level0.json");
        dicttype = Dict{String, Any},
    )
    V0, F0 = ESD._subdivide_icosahedron(Float64, 0)
    pairs = Tuple{Int,Int}[]
    for f in 1:size(F0, 2)
        a, b, c = F0[1, f], F0[2, f], F0[3, f]
        push!(pairs, REL.skolem_edge(a, b)); push!(pairs, REL.skolem_edge(b, c)); push!(pairs, REL.skolem_edge(c, a))
    end
    # Byte-for-byte canonical serialization (RFC §5.5.3) — the bytes every binding must match.
    @test REL.canonical_index_set_json(pairs) == golden["edge_set"]
    rk = REL.rank(pairs; base = 0)                       # canonical 0-based dense ids
    @test [rk.id[e] for e in rk.order] == Vector{Int}(golden["dense_ids_canonical"])
    @test length(rk.order) == golden["n_edges"]
end
