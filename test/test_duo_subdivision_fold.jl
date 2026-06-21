# DUO icosahedral build-time level-fold expressed declaratively (esd-heg.4, D3b).
#
# The full level-N mesh is the base icosahedron seed folded through the one-level
# value-invention refine pass (esd-heg.3 / subdivide_refine.esm) N times, where
# N = grid options.level is a build-time constant. This is the declarative
# replacement for the per-level driver loop `for _ in 1:level` of
# _subdivide_icosahedron (src/grids/duo.jl) — the last recursive kernel.
#
# The seed is the "literal/trivial FAQ constant": the 12 golden-ratio (phi)
# icosahedron vertices, normalized to the unit sphere by the same FAQ expression
# AST the refine pass uses (squares as products, 3-arg + left-assoc) so the seed
# is bit-for-bit identical to the imperative _icosahedron_vertices. The fold then
# chains the refine pass, feeding each pass's output vertices+faces as the next
# pass's input parameters.
#
# ESS exposes no native build-time repeat/iterate/fold construct over a build-time
# level (no such key in esm-schema.json; no whole-model value-invention runner), so
# the N-times unroll is performed caller-side here — exactly as the landed precedent
# test/test_duo_subdivision_faq.jl drives the single pass. The missing native
# affordance is filed as an ESS bead (see discretizations/grids/duo/faq/README.md).
#
# This test drives the *landed ESS engine primitives* (Relational.skolem_edge /
# rank for the value invention; eval_coeff for the seed + midpoint geometry FAQ —
# no shadow evaluator, per AGENTS.md single-pathway) and proves the fold reproduces
# the imperative _subdivide_icosahedron full mesh exactly for levels 0, 1, 2:
# imperative coords to ULP (bit-for-bit), faces position-for-position by coordinate
# (winding preserved), icos_level{0,1,2} counts (20*4^L cells), and the canonical
# full-mesh byte contract.

@testitem "duo build-time level-fold matches imperative full mesh" tags = [:grid, :duo, :faq, :subdivision, :fold] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const REL = EarthSciSerialization.Relational

    # --- Base icosahedron seed: the literal/trivial FAQ constant -------------
    # Raw golden-ratio vertices (pre-normalization). phi = (1+sqrt(5))/2.
    φ = (1.0 + sqrt(5.0)) / 2.0
    VERT_RAW = [
        0.0  1.0   φ;
        0.0 -1.0   φ;
        0.0  1.0  -φ;
        0.0 -1.0  -φ;
        1.0   φ  0.0;
       -1.0   φ  0.0;
        1.0  -φ  0.0;
       -1.0  -φ  0.0;
         φ  0.0  1.0;
         φ  0.0 -1.0;
        -φ  0.0  1.0;
        -φ  0.0 -1.0;
    ]
    FACE_VERT = [
        1  9  2; 1  2 11; 1 11  6; 1  6  5; 1  5  9; 9  5 10; 9 10  7; 9  7  2;
        2  7  8; 2  8 11; 11  8 12; 11 12  6; 6 12  3; 6  3  5; 5  3 10; 10  3  4;
        10  4  7; 7  4  8; 8  4 12; 12  4  3;
    ]

    # Seed normalize-to-sphere FAQ: coord_d = raw_d / sqrt(raw1² + raw2² + raw3²).
    # Squares are products (`*`), NOT `^`, and the 3-arg `+` folds left-assoc — the
    # exact float-op sequence of the imperative `n = sqrt(v[1]^2+v[2]^2+v[3]^2); v[d]/n`
    # (Julia `v^2` lowers to `v*v`). sqrt(Σ) is deterministic, so dividing each
    # component by it reproduces the single imperative `n` bit-for-bit.
    _rsq(n) = Dict{String,Any}("op" => "*", "args" => [n, n])
    _rnorm() = Dict{String,Any}("op" => "sqrt",
        "args" => [Dict{String,Any}("op" => "+", "args" => [_rsq("r1"), _rsq("r2"), _rsq("r3")])])
    _rcomp(d) = Dict{String,Any}("op" => "/", "args" => ["r$d", _rnorm()])

    function declarative_seed()
        Nv = size(VERT_RAW, 1)
        V = Matrix{Float64}(undef, 3, Nv)
        for i in 1:Nv
            b = Dict{String,Float64}("r1" => VERT_RAW[i, 1], "r2" => VERT_RAW[i, 2], "r3" => VERT_RAW[i, 3])
            V[1, i] = ESD.eval_coeff(_rcomp(1), b)
            V[2, i] = ESD.eval_coeff(_rcomp(2), b)
            V[3, i] = ESD.eval_coeff(_rcomp(3), b)
        end
        F = Matrix{Int}(permutedims(FACE_VERT))   # (3, 20)
        return V, F
    end

    # --- One refine pass (the esd-heg.3 value-invention FAQ) -----------------
    # Midpoint coordinate AST: (v_a[d] + v_b[d]) / sqrt(Σ_d (v_a[d]+v_b[d])²).
    _sum(d) = Dict{String,Any}("op" => "+", "args" => ["a$d", "b$d"])
    _sq(d) = Dict{String,Any}("op" => "*", "args" => [_sum(d), _sum(d)])
    _norm() = Dict{String,Any}("op" => "sqrt",
        "args" => [Dict{String,Any}("op" => "+", "args" => [_sq(1), _sq(2), _sq(3)])])
    _component(d) = Dict{String,Any}("op" => "/", "args" => [_sum(d), _norm()])

    function faq_refine_once(V::Matrix{Float64}, F::Matrix{Int})
        Nv = size(V, 2)
        Nc = size(F, 2)
        pairs = Tuple{Int,Int}[]
        for f in 1:Nc
            a, b, c = F[1, f], F[2, f], F[3, f]
            push!(pairs, REL.skolem_edge(a, b))
            push!(pairs, REL.skolem_edge(b, c))
            push!(pairs, REL.skolem_edge(c, a))
        end
        rk = REL.rank(pairs; base = 1)          # rk.order = canonical sorted edges; rk.id = dense id
        edges = rk.order
        Ne = length(edges)

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

        midvid(x, y) = Nv + rk.id[REL.skolem_edge(x, y)]

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

    # The build-time level-fold: seed, then chain the refine pass `level` times.
    function fold_from_seed(level::Int)
        V, F = declarative_seed()
        for _ in 1:level
            V, F, _ = faq_refine_once(V, F)
        end
        return V, F
    end

    coordset(V) = Set(Tuple(V[:, j]) for j in 1:size(V, 2))
    facetris(V, F) = [(Tuple(V[:, F[1, k]]), Tuple(V[:, F[2, k]]), Tuple(V[:, F[3, k]])) for k in 1:size(F, 2)]

    # (A) The declarative φ-seed reproduces the imperative seed bit-for-bit — the
    #     fold starts from the identical level-0 mesh, not from _subdivide_icosahedron.
    Vseed, Fseed = declarative_seed()
    V0imp, F0imp = ESD._subdivide_icosahedron(Float64, 0)
    @test Vseed == V0imp
    @test Fseed == F0imp

    # (B) Folding the refine pass `level` times reproduces the imperative full mesh,
    #     levels 0, 1, 2.
    for L in 0:2
        Vf, Ff = fold_from_seed(L)
        Vexp, Fexp = ESD._subdivide_icosahedron(Float64, L)

        # Counts: 20·4^L cells, 10·4^L+2 vertices, 30·4^L edges.
        @test size(Vf, 2) == size(Vexp, 2) == 10 * 4^L + 2
        @test size(Ff, 2) == size(Fexp, 2) == 20 * 4^L

        # Same vertex coordinate set, bit-for-bit (seed FAQ + midpoint FAQ ≡ imperative).
        @test coordset(Vf) == coordset(Vexp)

        # Same triangles, position-for-position with winding preserved — only the
        # invented midpoint *labels* differ (§5.7 canonical vs imperative push! order).
        @test facetris(Vf, Ff) == facetris(Vexp, Fexp)
    end

    # (C) Determinism of the value-invented labeling: the §5.7 canonical vertex
    #     order (base verts, then midpoints in distinct+rank order) is independent of
    #     seed face iteration order AND triangle winding — `skolem_edge` is undirected,
    #     so the invented edge set at every level is identical regardless of how the
    #     faces are presented. Hence the folded *vertex array* is bit-for-bit identical.
    #     (Determinism is asserted on the invented vertices, as in esd-heg.3's edge-set
    #     test: the face triangles are winding-dependent by construction — reversing a
    #     seed face's winding reverses its children — so the face set is NOT, and need
    #     not be, winding-invariant.)
    function fold_verts(V, F, level)
        for _ in 1:level
            V, F, _ = faq_refine_once(V, F)
        end
        return V
    end
    Fperm = Fseed[:, reverse(1:size(Fseed, 2))]                     # permuted face order
    Frev = vcat(Fseed[3:3, :], Fseed[2:2, :], Fseed[1:1, :])        # reversed winding
    Vf_norm = fold_verts(Vseed, Fseed, 2)
    @test fold_verts(copy(Vseed), Fperm, 2) == Vf_norm
    @test fold_verts(copy(Vseed), Frev, 2) == Vf_norm
end

@testitem "duo level-fold: declarative seed doc + canonical full-mesh byte contract" tags = [:grid, :duo, :faq, :subdivision, :fold, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const REL = EarthSciSerialization.Relational
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    DUO_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo")
    FAQ_DIR = joinpath(DUO_DIR, "faq")

    # --- The declarative seed/fold document exists, is schema-valid, and expresses
    #     the seed normalize-to-sphere FAQ.
    doc = JSON.parsefile(joinpath(FAQ_DIR, "subdivide_fold.esm"); dicttype = Dict{String, Any})
    model = doc["models"]["IcosahedralSeed"]
    isets = model["index_sets"]
    @test isets["verts"]["size"] == 12
    @test isets["faces"]["size"] == 20
    vars = model["variables"]
    @test vars["vert_raw"]["type"] == "parameter"
    @test vars["face_vert"]["type"] == "parameter"
    @test vars["seed_coord"]["type"] == "state"
    eqs = model["equations"]
    # seed_coord = vert_raw / sqrt(seed_sq): a normalize-to-sphere FAQ (division by sqrt).
    seed_coord_eq = only(filter(e -> e["lhs"]["args"][1] == "seed_coord", eqs))
    @test seed_coord_eq["rhs"]["expr"]["op"] == "/"
    @test any(e -> e["lhs"]["args"][1] == "seed_sq", eqs)
    # Schema validity (mirrors subdivide_refine.esm; the claim the README makes).
    @test isempty(EarthSciSerialization.validate_schema(doc))

    # --- Rebuild the fold from the declarative seed, drive the landed ESS primitives.
    φ = (1.0 + sqrt(5.0)) / 2.0
    VERT_RAW = [
        0.0  1.0   φ;  0.0 -1.0   φ;  0.0  1.0  -φ;  0.0 -1.0  -φ;
        1.0   φ  0.0; -1.0   φ  0.0;  1.0  -φ  0.0; -1.0  -φ  0.0;
         φ  0.0  1.0;   φ  0.0 -1.0;  -φ  0.0  1.0;  -φ  0.0 -1.0;
    ]
    FACE_VERT = [
        1  9  2; 1  2 11; 1 11  6; 1  6  5; 1  5  9; 9  5 10; 9 10  7; 9  7  2;
        2  7  8; 2  8 11; 11  8 12; 11 12  6; 6 12  3; 6  3  5; 5  3 10; 10  3  4;
        10  4  7; 7  4  8; 8  4 12; 12  4  3;
    ]
    _rsq(n) = Dict{String,Any}("op" => "*", "args" => [n, n])
    _rnorm() = Dict{String,Any}("op" => "sqrt",
        "args" => [Dict{String,Any}("op" => "+", "args" => [_rsq("r1"), _rsq("r2"), _rsq("r3")])])
    _rcomp(d) = Dict{String,Any}("op" => "/", "args" => ["r$d", _rnorm()])
    _sum(d) = Dict{String,Any}("op" => "+", "args" => ["a$d", "b$d"])
    _sq(d) = Dict{String,Any}("op" => "*", "args" => [_sum(d), _sum(d)])
    _norm() = Dict{String,Any}("op" => "sqrt",
        "args" => [Dict{String,Any}("op" => "+", "args" => [_sq(1), _sq(2), _sq(3)])])
    _component(d) = Dict{String,Any}("op" => "/", "args" => [_sum(d), _norm()])

    function seed_mesh()
        V = Matrix{Float64}(undef, 3, size(VERT_RAW, 1))
        for i in 1:size(VERT_RAW, 1)
            b = Dict{String,Float64}("r1" => VERT_RAW[i, 1], "r2" => VERT_RAW[i, 2], "r3" => VERT_RAW[i, 3])
            V[1, i] = ESD.eval_coeff(_rcomp(1), b); V[2, i] = ESD.eval_coeff(_rcomp(2), b); V[3, i] = ESD.eval_coeff(_rcomp(3), b)
        end
        return V, Matrix{Int}(permutedims(FACE_VERT))
    end
    function refine_once(V, F)
        Nv = size(V, 2); Nc = size(F, 2)
        pairs = Tuple{Int,Int}[]
        for f in 1:Nc
            a, b, c = F[1, f], F[2, f], F[3, f]
            push!(pairs, REL.skolem_edge(a, b)); push!(pairs, REL.skolem_edge(b, c)); push!(pairs, REL.skolem_edge(c, a))
        end
        rk = REL.rank(pairs; base = 1); edges = rk.order; Ne = length(edges)
        Vnew = Matrix{Float64}(undef, 3, Nv + Ne); Vnew[:, 1:Nv] = V
        for (i, (lo, hi)) in enumerate(edges)
            b = Dict{String,Float64}("a1" => V[1, lo], "a2" => V[2, lo], "a3" => V[3, lo],
                "b1" => V[1, hi], "b2" => V[2, hi], "b3" => V[3, hi])
            Vnew[1, Nv + i] = ESD.eval_coeff(_component(1), b); Vnew[2, Nv + i] = ESD.eval_coeff(_component(2), b); Vnew[3, Nv + i] = ESD.eval_coeff(_component(3), b)
        end
        midvid(x, y) = Nv + rk.id[REL.skolem_edge(x, y)]
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
        return Vnew, Fnew
    end
    function fold(level)
        V, F = seed_mesh()
        for _ in 1:level; V, F = refine_once(V, F); end
        return V, F
    end
    faces_json(F) = "[" * join(["[$(F[1,k]),$(F[2,k]),$(F[3,k])]" for k in 1:size(F, 2)], ",") * "]"

    for L in 0:2
        _, Ff = fold(L)

        # Counts come from the committed icos_level{L}.esm grid descriptors.
        desc = JSON.parsefile(joinpath(DUO_DIR, "icos_level$(L).esm"); dicttype = Dict{String, Any})
        @test size(Ff, 2) == desc["n_cells"] == 20 * 4^L

        # Canonical full-mesh byte contract: the §5.7-labeled faces every binding
        # folding from the canonical seed must reproduce.
        golden = JSON.parsefile(
            joinpath(FAQ_DIR, "fixtures", "canonical", "mesh_level$(L).json");
            dicttype = Dict{String, Any},
        )
        @test golden["n_cells"] == desc["n_cells"]
        @test golden["n_vertices"] == desc["n_vertices"]
        @test golden["n_edges"] == desc["n_edges"]
        @test faces_json(Ff) == golden["faces_canonical"]
    end
end
