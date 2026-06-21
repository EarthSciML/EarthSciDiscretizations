# DUO edge + MPAS Voronoi-dual GEOMETRY expressed as a value-invention FAQ (esd-heg.7 / D2b).
#
# The geometry half of the imperative `_duo_voronoi_dual` (src/grids/mpas.jl, steps
# 1 / circumcenters / 5 / area) and the DUO Tier-U `edge_length` / `cell_distance`
# (src/grids/duo.jl) is the RFC `semiring-faq-unified-ir` §8.1 geometry-FAQ pattern:
# elementwise circumcenters (the cross-product sum a×b+b×c+c×a, normalized),
# per-edge great-circle arcs `R·acos(·)`, the geographic `atan2`/`asin` of cell
# centres and edge midpoints, and the Voronoi dual-cell area as a spherical-excess
# (L'Huilier) FAN over the angular circumcenter ring. The declarative spec is
# `discretizations/grids/duo/faq/edge_dual_geometry.esm`; this test drives the
# *landed ESS engine* (`eval_coeff` — the single-pathway passthrough, no shadow
# evaluator per AGENTS.md) on the same expressions and proves they reproduce the
# imperative geometry *bit-for-bit* (strict Float64 bit identity) across the dual
# level ladder.
#
# Why bit-exact and not merely close:
#   * Inputs are the imperative's recovered-unit coordinates `duo.vertices ./ R` and
#     `duo.cell_cart ./ R` (not the raw subdivision unit verts), so the float ops are
#     fed identical operands.
#   * Squares are products (`x*x`, not `^`); dot products and the centroid/midpoint
#     sums fold in space order; the L'Huilier tan-product is left-folded.
#   * The imperative circumcenter applies a universal outward-orientation flip
#     (`cc·centroid < 0 → negate`, which fires for EVERY face on the consistently
#     wound icosahedral mesh). That is reproduced as the standard cross sum NEGATED
#     (×−1.0) — bit-for-bit, INCLUDING the `−0.0` a negated zero component carries
#     (reversing the operand order would give `x−x = +0.0` instead and break the
#     contract on the symmetric faces).
#   * The only divergence from the imperative is dropping the `clamp`/`max` clipping
#     guards (acceptance: deterministic formulas, no clipping): on a valid icosahedral
#     mesh every `acos` argument and the L'Huilier radicand already lie in range, so
#     clipping never fires and the bytes match — confirmed at levels 1–3.
#
# The angularly-ordered incident-face ring of each dual cell is the D1b
# `voronoi_dual_topology_faq(...).sorted_vertex_faces` (esd-heg.2), exposed exactly
# so this geometry FAQ can fan the dual-cell area over it.

@testitem "duo edge+dual geometry via FAQ matches imperative (bit-identical)" tags = [:grid, :duo, :mpas, :faq, :geometry] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization

    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    rawbits(x) = reinterpret(Int64, Float64(x))
    same(got, ref) = rawbits(got) == rawbits(ref)   # strict bitwise (catches ±0.0)

    # --- circumcenter direction = −1.0 · (standard cross sum a×b+b×c+c×a):
    #     a=corner1, b=corner2, c=corner3; component k: 1=x,2=y,3=z.
    negc(t1, t2, t3) = mk("*", -1.0, mk("+", t1, t2, t3))
    ccx = negc(mk("-", mk("*", "a2", "b3"), mk("*", "a3", "b2")),
        mk("-", mk("*", "b2", "c3"), mk("*", "b3", "c2")),
        mk("-", mk("*", "c2", "a3"), mk("*", "c3", "a2")))
    ccy = negc(mk("-", mk("*", "a3", "b1"), mk("*", "a1", "b3")),
        mk("-", mk("*", "b3", "c1"), mk("*", "b1", "c3")),
        mk("-", mk("*", "c3", "a1"), mk("*", "c1", "a3")))
    ccz = negc(mk("-", mk("*", "a1", "b2"), mk("*", "a2", "b1")),
        mk("-", mk("*", "b1", "c2"), mk("*", "b2", "c1")),
        mk("-", mk("*", "c1", "a2"), mk("*", "c2", "a1")))

    # --- great-circle arc R·acos(p·q) over point leaves p1..p3, q1..q3 (no clamp).
    arc = mk("*", "R", mk("acos",
        mk("+", mk("*", "p1", "q1"), mk("*", "p2", "q2"), mk("*", "p3", "q3"))))

    # --- spherical-excess (L'Huilier) of triangle (A,B,C), NO R² (mirrors
    #     _spherical_triangle_area): da=acos(B·C), db=acos(C·A), dc=acos(A·B).
    dot_ab = mk("+", mk("*", "a1", "b1"), mk("*", "a2", "b2"), mk("*", "a3", "b3"))
    dot_bc = mk("+", mk("*", "b1", "c1"), mk("*", "b2", "c2"), mk("*", "b3", "c3"))
    dot_ca = mk("+", mk("*", "c1", "a1"), mk("*", "c2", "a2"), mk("*", "c3", "a3"))
    sda = mk("acos", dot_bc); sdb = mk("acos", dot_ca); sdc = mk("acos", dot_ab)
    ss = mk("*", 0.5, mk("+", sda, sdb, sdc))
    shalf(x) = mk("/", x, 2.0)
    st = mk("*", mk("tan", shalf(ss)), mk("tan", shalf(mk("-", ss, sda))),
        mk("tan", shalf(mk("-", ss, sdb))), mk("tan", shalf(mk("-", ss, sdc))))
    excess = mk("*", 4.0, mk("atan", mk("sqrt", st)))

    # --- lon/lat of a unit point u1..u3, and of a normalized edge midpoint.
    lon_u = mk("atan2", "u2", "u1"); lat_u = mk("asin", "u3")
    mxx = mk("+", "a1", "b1"); myy = mk("+", "a2", "b2"); mzz = mk("+", "a3", "b3")
    mnn = mk("sqrt", mk("+", mk("*", mxx, mxx), mk("*", myy, myy), mk("*", mzz, mzz)))
    lon_e = mk("atan2", mk("/", myy, mnn), mk("/", mxx, mnn))
    lat_e = mk("asin", mk("/", mzz, mnn))

    cbind(V, v1, v2, v3) = Dict{String, Float64}(
        "a1" => V[1, v1], "a2" => V[2, v1], "a3" => V[3, v1],
        "b1" => V[1, v2], "b2" => V[2, v2], "b3" => V[3, v2],
        "c1" => V[1, v3], "c2" => V[2, v3], "c3" => V[3, v3])
    pq(P, i, Q, j, R) = Dict{String, Float64}(
        "p1" => P[1, i], "p2" => P[2, i], "p3" => P[3, i],
        "q1" => Q[1, j], "q2" => Q[2, j], "q3" => Q[3, j], "R" => R)

    R = 6.371e6
    # The MPAS Voronoi dual requires a subdivided primal (level ≥ 1).
    for level in 1:3
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        imp = ESD._duo_voronoi_dual(level; R = R)
        topo = voronoi_dual_topology_faq(g.vertices, g.faces, g.cell_cart; R = R)
        svf = topo.sorted_vertex_faces
        ec = ESD._duo_edge_cells(g)             # edge → the two incident DUO faces
        Nv = size(g.vertices, 2); Nc = size(g.faces, 2); Ne = size(g.edges, 2)
        Rf = Float64(g.R)
        vu = Float64.(g.vertices) ./ Rf         # recovered-unit cell centres (DUO verts)
        cu = Float64.(g.cell_cart) ./ Rf        # recovered-unit face centroids
        @test Nv == imp.n_cells

        nbad = 0

        # (1) circumcenters: FAQ (standard cross sum × −1.0) vs imperative face_cc.
        facecc = Matrix{Float64}(undef, 3, Nc)
        for f in 1:Nc
            v1, v2, v3 = g.faces[1, f], g.faces[2, f], g.faces[3, f]
            a1x, a1y, a1z = vu[1, v1], vu[2, v1], vu[3, v1]
            a2x, a2y, a2z = vu[1, v2], vu[2, v2], vu[3, v2]
            a3x, a3y, a3z = vu[1, v3], vu[2, v3], vu[3, v3]
            ix = (a1y * a2z - a1z * a2y) + (a2y * a3z - a2z * a3y) + (a3y * a1z - a3z * a1y)
            iy = (a1z * a2x - a1x * a2z) + (a2z * a3x - a2x * a3z) + (a3z * a1x - a3x * a1z)
            iz = (a1x * a2y - a1y * a2x) + (a2x * a3y - a2y * a3x) + (a3x * a1y - a3y * a1x)
            if ix * cu[1, f] + iy * cu[2, f] + iz * cu[3, f] < 0   # the universal flip
                ix = -ix; iy = -iy; iz = -iz
            end
            cn = sqrt(ix^2 + iy^2 + iz^2)
            iccx, iccy, iccz = ix / cn, iy / cn, iz / cn
            b = cbind(vu, v1, v2, v3)
            rx = ESD.eval_coeff(ccx, b); ry = ESD.eval_coeff(ccy, b); rz = ESD.eval_coeff(ccz, b)
            fn = sqrt(rx * rx + ry * ry + rz * rz)
            facecc[1, f] = rx / fn; facecc[2, f] = ry / fn; facecc[3, f] = rz / fn
            nbad += !same(facecc[1, f], iccx) + !same(facecc[2, f], iccy) + !same(facecc[3, f], iccz)
        end

        # (2) dc_edge (cell↔cell arc) — also the DUO Tier-U edge_length.
        el = edge_length(g)
        for e in 1:Ne
            v1, v2 = g.edges[1, e], g.edges[2, e]
            v = ESD.eval_coeff(arc, pq(vu, v1, vu, v2, Rf))
            nbad += !same(v, imp.dc_edge[e]) + !same(v, el[e])
        end

        # (3) dv_edge (circumcenter↔circumcenter arc).
        for e in 1:Ne
            f1, f2 = ec[1, e], ec[2, e]
            v = ESD.eval_coeff(arc, pq(facecc, f1, facecc, f2, Rf))
            nbad += !same(v, imp.dv_edge[e])
        end

        # (4) cell_distance (DUO Tier-U centroid↔centroid arc).
        cd = cell_distance(g)
        for e in 1:Ne
            f1, f2 = ec[1, e], ec[2, e]
            v = ESD.eval_coeff(arc, pq(cu, f1, cu, f2, Rf))
            nbad += !same(v, cd[e])
        end

        # (5) area_cell — spherical-excess fan over the angular circumcenter ring.
        for v in 1:Nv
            sfaces = svf[v]; kv = length(sfaces)
            area_v = 0.0
            for i in 1:kv
                fi = sfaces[i]; fnext = sfaces[(i % kv) + 1]
                b = Dict{String, Float64}(
                    "a1" => vu[1, v], "a2" => vu[2, v], "a3" => vu[3, v],
                    "b1" => facecc[1, fi], "b2" => facecc[2, fi], "b3" => facecc[3, fi],
                    "c1" => facecc[1, fnext], "c2" => facecc[2, fnext], "c3" => facecc[3, fnext])
                area_v += ESD.eval_coeff(excess, b)
            end
            nbad += !same(Rf * Rf * area_v, imp.area_cell[v])
        end

        # (6) lon/lat of cells and edges.
        for v in 1:Nv
            b = Dict{String, Float64}("u1" => vu[1, v], "u2" => vu[2, v], "u3" => vu[3, v])
            nbad += !same(ESD.eval_coeff(lon_u, b), imp.lon_cell[v]) +
                !same(ESD.eval_coeff(lat_u, b), imp.lat_cell[v])
        end
        for e in 1:Ne
            v1, v2 = g.edges[1, e], g.edges[2, e]
            b = Dict{String, Float64}(
                "a1" => vu[1, v1], "a2" => vu[2, v1], "a3" => vu[3, v1],
                "b1" => vu[1, v2], "b2" => vu[2, v2], "b3" => vu[3, v2])
            nbad += !same(ESD.eval_coeff(lon_e, b), imp.lon_edge[e]) +
                !same(ESD.eval_coeff(lat_e, b), imp.lat_edge[e])
        end

        @test nbad == 0   # bit-for-bit, every face/edge/cell, every metric
    end

    # Physical sanity: the FAQ dual-cell areas tile the full sphere 4πR² and each is
    # positive — independent confirmation the geometry is right, not just self-consistent.
    let level = 1
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        topo = voronoi_dual_topology_faq(g.vertices, g.faces, g.cell_cart; R = R)
        svf = topo.sorted_vertex_faces
        Nc = size(g.faces, 2); Nv = size(g.vertices, 2)
        Rf = Float64(g.R)
        vu = Float64.(g.vertices) ./ Rf; cu = Float64.(g.cell_cart) ./ Rf
        facecc = Matrix{Float64}(undef, 3, Nc)
        for f in 1:Nc
            v1, v2, v3 = g.faces[1, f], g.faces[2, f], g.faces[3, f]
            b = cbind(vu, v1, v2, v3)
            rx = ESD.eval_coeff(ccx, b); ry = ESD.eval_coeff(ccy, b); rz = ESD.eval_coeff(ccz, b)
            fn = sqrt(rx * rx + ry * ry + rz * rz)
            facecc[1, f] = rx / fn; facecc[2, f] = ry / fn; facecc[3, f] = rz / fn
        end
        faq_area = Float64[]
        for v in 1:Nv
            sfaces = svf[v]; kv = length(sfaces)
            area_v = 0.0
            for i in 1:kv
                fi = sfaces[i]; fnext = sfaces[(i % kv) + 1]
                b = Dict{String, Float64}(
                    "a1" => vu[1, v], "a2" => vu[2, v], "a3" => vu[3, v],
                    "b1" => facecc[1, fi], "b2" => facecc[2, fi], "b3" => facecc[3, fi],
                    "c1" => facecc[1, fnext], "c2" => facecc[2, fnext], "c3" => facecc[3, fnext])
                area_v += ESD.eval_coeff(excess, b)
            end
            push!(faq_area, Rf * Rf * area_v)
        end
        @test all(>(0), faq_area)
        @test sum(faq_area) ≈ 4π * Rf^2 rtol = 1e-12
    end
end

@testitem "duo edge+dual geometry FAQ: declarative doc + cross-binding byte contract" tags = [:grid, :duo, :mpas, :faq, :geometry, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const ESS = EarthSciSerialization
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    FAQ_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo", "faq")

    # --- The declarative value-invention FAQ document exists, is schema-valid, and
    #     expresses the dual-geometry chain.
    doc_path = joinpath(FAQ_DIR, "edge_dual_geometry.esm")
    doc = JSON.parsefile(doc_path; dicttype = Dict{String, Any})
    @test isempty(ESS.validate_schema(doc))   # validates against the bundled esm-schema.json

    model = doc["models"]["DuoEdgeDualGeometry"]
    isets = model["index_sets"]
    for s in ("cells", "faces", "edges", "ring", "corners", "space", "ends")
        @test haskey(isets, s)
    end
    @test isets["corners"]["size"] == 3 && isets["space"]["size"] == 3 && isets["ends"]["size"] == 2

    vars = model["variables"]
    for v in ("R", "vert_coord", "centroid", "face_vert", "edge_cell", "edge_face",
            "ring_face", "ring_face_next", "circ_x", "circ_y", "circ_z",
            "dc_edge", "dv_edge", "cell_distance", "area_cell",
            "lon_cell", "lat_cell", "lon_edge", "lat_edge")
        @test haskey(vars, v)
    end
    for p in ("R", "vert_coord", "centroid", "face_vert", "edge_cell", "edge_face", "ring_face", "ring_face_next")
        @test vars[p]["type"] == "parameter"
    end
    for s in ("circ_x", "dc_edge", "dv_edge", "cell_distance", "area_cell", "lon_cell", "lon_edge")
        @test vars[s]["type"] == "state"
    end

    eqs = model["equations"]
    lhs_name(e) = (e["lhs"] isa Dict && get(e["lhs"], "op", "") == "index") ? e["lhs"]["args"][1] : nothing
    eqfor(name) = only(filter(e -> lhs_name(e) == name, eqs))["rhs"]
    # Every named output is produced by an aggregate (the geometry FAQ).
    for name in ("circ_x", "dc_edge", "dv_edge", "cell_distance", "area_cell",
            "lon_cell", "lat_cell", "lon_edge", "lat_edge")
        @test eqfor(name)["op"] == "aggregate"
    end
    # The dual-cell area reduces over the ring (the spherical-excess fan).
    @test haskey(eqfor("area_cell")["ranges"], "s")
    @test eqfor("area_cell")["ranges"]["s"]["from"] == "ring"

    # Recursively collect operator names used by a node's rhs (ignoring prose).
    function ops_of(node, acc = Set{String}())
        if node isa AbstractDict
            haskey(node, "op") && push!(acc, String(node["op"]))
            for (k, v) in node
                k == "_comment" || ops_of(v, acc)
            end
        elseif node isa AbstractVector
            for v in node
                ops_of(v, acc)
            end
        end
        return acc
    end
    # Arcs are R·acos(dot); lon/lat use atan2/asin; the area fan uses the spherical-excess
    # leaves (tan/atan/sqrt); the circumcenter cross product folds via +/−/*.
    @test "acos" in ops_of(eqfor("dc_edge"))
    @test "acos" in ops_of(eqfor("dv_edge"))
    @test "acos" in ops_of(eqfor("cell_distance"))
    @test "atan2" in ops_of(eqfor("lon_cell"))
    @test "asin" in ops_of(eqfor("lat_cell"))
    @test "atan2" in ops_of(eqfor("lon_edge"))
    let a = ops_of(eqfor("fan_excess"))
        @test "tan" in a && "atan" in a && "sqrt" in a
    end
    let cc = ops_of(eqfor("ccr_x"))
        @test "+" in cc && "-" in cc && "*" in cc       # cross-product sum (negated)
    end
    # Deterministic formulas, no clipping (acceptance contract): no clamp/clip/max OP in
    # any equation, and squares are products (no `^`).
    all_ops = reduce((a, e) -> ops_of(e["rhs"], a), eqs; init = Set{String}())
    @test isempty(intersect(all_ops, Set(["clamp", "clip", "max", "^"])))

    # --- Cross-binding byte/value contract: the pinned level-1 golden is reproduced
    #     bit-for-bit by BOTH the imperative builder and the FAQ (ESS eval_coeff).
    golden = JSON.parsefile(
        joinpath(FAQ_DIR, "fixtures", "canonical", "edge_dual_geometry_level1.json");
        dicttype = Dict{String, Any},
    )
    R = 6.371e6
    level = golden["level"]
    g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
    imp = ESD._duo_voronoi_dual(level; R = R)
    topo = voronoi_dual_topology_faq(g.vertices, g.faces, g.cell_cart; R = R)
    ec = ESD._duo_edge_cells(g)
    Nv = size(g.vertices, 2); Nc = size(g.faces, 2); Ne = size(g.edges, 2)
    Rf = Float64(g.R)
    vu = Float64.(g.vertices) ./ Rf; cu = Float64.(g.cell_cart) ./ Rf
    @test golden["n_cells"] == Nv && golden["n_faces"] == Nc && golden["n_edges"] == Ne
    rb(x) = reinterpret(Int64, Float64(x))

    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    negc(t1, t2, t3) = mk("*", -1.0, mk("+", t1, t2, t3))
    ccx = negc(mk("-", mk("*", "a2", "b3"), mk("*", "a3", "b2")),
        mk("-", mk("*", "b2", "c3"), mk("*", "b3", "c2")),
        mk("-", mk("*", "c2", "a3"), mk("*", "c3", "a2")))
    ccy = negc(mk("-", mk("*", "a3", "b1"), mk("*", "a1", "b3")),
        mk("-", mk("*", "b3", "c1"), mk("*", "b1", "c3")),
        mk("-", mk("*", "c3", "a1"), mk("*", "c1", "a3")))
    ccz = negc(mk("-", mk("*", "a1", "b2"), mk("*", "a2", "b1")),
        mk("-", mk("*", "b1", "c2"), mk("*", "b2", "c1")),
        mk("-", mk("*", "c1", "a2"), mk("*", "c2", "a1")))
    arc = mk("*", "R", mk("acos",
        mk("+", mk("*", "p1", "q1"), mk("*", "p2", "q2"), mk("*", "p3", "q3"))))
    cbind(V, v1, v2, v3) = Dict{String, Float64}(
        "a1" => V[1, v1], "a2" => V[2, v1], "a3" => V[3, v1],
        "b1" => V[1, v2], "b2" => V[2, v2], "b3" => V[3, v2],
        "c1" => V[1, v3], "c2" => V[2, v3], "c3" => V[3, v3])
    pq(P, i, Q, j) = Dict{String, Float64}(
        "p1" => P[1, i], "p2" => P[2, i], "p3" => P[3, i],
        "q1" => Q[1, j], "q2" => Q[2, j], "q3" => Q[3, j], "R" => Rf)

    # FAQ circumcenters (standard cross sum × −1.0) — and imperative reference.
    facecc = Matrix{Float64}(undef, 3, Nc)
    for f in 1:Nc
        v1, v2, v3 = g.faces[1, f], g.faces[2, f], g.faces[3, f]
        a1x, a1y, a1z = vu[1, v1], vu[2, v1], vu[3, v1]
        a2x, a2y, a2z = vu[1, v2], vu[2, v2], vu[3, v2]
        a3x, a3y, a3z = vu[1, v3], vu[2, v3], vu[3, v3]
        ix = (a1y * a2z - a1z * a2y) + (a2y * a3z - a2z * a3y) + (a3y * a1z - a3z * a1y)
        iy = (a1z * a2x - a1x * a2z) + (a2z * a3x - a2x * a3z) + (a3z * a1x - a3x * a1z)
        iz = (a1x * a2y - a1y * a2x) + (a2x * a3y - a2y * a3x) + (a3x * a1y - a3y * a1x)
        if ix * cu[1, f] + iy * cu[2, f] + iz * cu[3, f] < 0
            ix = -ix; iy = -iy; iz = -iz
        end
        icn = sqrt(ix^2 + iy^2 + iz^2)
        impx, impy, impz = ix / icn, iy / icn, iz / icn
        b = cbind(vu, v1, v2, v3)
        rx = ESD.eval_coeff(ccx, b); ry = ESD.eval_coeff(ccy, b); rz = ESD.eval_coeff(ccz, b)
        fn = sqrt(rx * rx + ry * ry + rz * rz)
        facecc[1, f] = rx / fn; facecc[2, f] = ry / fn; facecc[3, f] = rz / fn
        # (a) imperative == golden; (b) FAQ == golden — bit-for-bit, both paths.
        @test rb(golden["circ_x"][f]) == rb(impx) == rb(facecc[1, f])
        @test rb(golden["circ_y"][f]) == rb(impy) == rb(facecc[2, f])
        @test rb(golden["circ_z"][f]) == rb(impz) == rb(facecc[3, f])
    end

    el = edge_length(g); cd = cell_distance(g)
    for e in 1:Ne
        v1, v2 = g.edges[1, e], g.edges[2, e]
        f1, f2 = ec[1, e], ec[2, e]
        faq_dc = ESD.eval_coeff(arc, pq(vu, v1, vu, v2))
        faq_dv = ESD.eval_coeff(arc, pq(facecc, f1, facecc, f2))
        faq_cd = ESD.eval_coeff(arc, pq(cu, f1, cu, f2))
        @test rb(golden["dc_edge"][e]) == rb(imp.dc_edge[e]) == rb(faq_dc)
        @test rb(golden["dv_edge"][e]) == rb(imp.dv_edge[e]) == rb(faq_dv)
        @test rb(golden["cell_distance"][e]) == rb(cd[e]) == rb(faq_cd)
        @test rb(el[e]) == rb(imp.dc_edge[e])   # edge_length ≡ dc_edge
    end

    for v in 1:Nv
        @test rb(golden["lon_cell"][v]) == rb(imp.lon_cell[v])
        @test rb(golden["lat_cell"][v]) == rb(imp.lat_cell[v])
        @test rb(golden["area_cell"][v]) == rb(imp.area_cell[v])
    end
    for e in 1:Ne
        @test rb(golden["lon_edge"][e]) == rb(imp.lon_edge[e])
        @test rb(golden["lat_edge"][e]) == rb(imp.lat_edge[e])
    end
end
