# DUO icosahedral primal-cell geometry expressed as a value-invention FAQ (esd-heg.6).
#
# The per-cell geometry loop in `build_duo_grid` (src/grids/duo.jl) — cell area
# (spherical excess, L'Huilier), the cell-center centroid (normalized mean of the
# three corner unit vectors), and the geographic lon/lat (atan2/asin) — is the
# RFC `semiring-faq-unified-ir` §8.1 geometry-FAQ pattern: per cell, a `sum_product`
# aggregate over the three corner vertices producing each scalar metric from the
# scalar leaves `acos`/`asin`/`atan`/`atan2`/`sqrt`/`tan` (ESS bead ess-9x1, the P0
# inverse-trig leaves). The declarative spec is
# `discretizations/grids/duo/faq/primal_geometry.esm`; this test drives the *landed
# ESS engine* (`eval_coeff` — the single-pathway passthrough, no shadow evaluator per
# AGENTS.md) on the same expressions and proves they reproduce the imperative geometry
# *bit-for-bit* (0 ULP) across the level ladder.
#
# Why bit-exact and not merely close: the FAQ mirrors the imperative float ops —
# squares are products (`x*x`, not `^`); the centroid sum, the squared norm, and the
# three side-arc dot products fold in corner/space order; the L'Huilier tan-product is
# left-folded. The only intentional divergence from the imperative code is the removal
# of the `clamp`/`max` clipping guards (the acceptance asks for deterministic formulas,
# no clipping): on a valid icosahedral mesh every `acos`/`asin` argument and the
# L'Huilier radicand already lie in range, so clipping never fires and the bytes match.

@testitem "duo primal geometry via FAQ matches imperative (0 ULP)" tags = [:grid, :duo, :faq, :geometry] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization

    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

    # --- The geometry ASTs over scalar corner leaves a1..a3 (corner 1), b1..b3
    #     (corner 2), c1..c3 (corner 3) and the radius R. These are the scalar
    #     instances of the indexed aggregates in primal_geometry.esm, and they mirror
    #     build_duo_grid's per-cell loop exactly (minus clipping).

    # centroid: m = a+b+c per component; n = sqrt(mx^2+my^2+mz^2); u = m/n.
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nsq = mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz))   # squares as *, not ^
    nrm = mk("sqrt", nsq)
    ux = mk("/", mx, nrm); uy = mk("/", my, nrm); uz = mk("/", mz, nrm)
    cart = (mk("*", "R", ux), mk("*", "R", uy), mk("*", "R", uz))
    lon = mk("atan2", uy, ux)
    lat = mk("asin", uz)                                              # no clamp

    # area (L'Huilier): da=acos(b·c), db=acos(c·a), dc=acos(a·b); s=(da+db+dc)/2;
    # t=tan(s/2)·tan((s-da)/2)·tan((s-db)/2)·tan((s-dc)/2); area = 4·atan(sqrt t)·R².
    dot_bc = mk("+", mk("*", "b1", "c1"), mk("*", "b2", "c2"), mk("*", "b3", "c3"))
    dot_ca = mk("+", mk("*", "c1", "a1"), mk("*", "c2", "a2"), mk("*", "c3", "a3"))
    dot_ab = mk("+", mk("*", "a1", "b1"), mk("*", "a2", "b2"), mk("*", "a3", "b3"))
    da = mk("acos", dot_bc); db = mk("acos", dot_ca); dc = mk("acos", dot_ab)
    s = mk("*", 0.5, mk("+", da, db, dc))
    half(x) = mk("/", x, 2.0)
    t = mk("*", mk("tan", half(s)), mk("tan", half(mk("-", s, da))),
        mk("tan", half(mk("-", s, db))), mk("tan", half(mk("-", s, dc))))
    area = mk("*", mk("*", 4.0, mk("atan", mk("sqrt", t))), mk("*", "R", "R"))   # no max guard

    bindings(V, i, j, k, R) = Dict{String, Float64}(
        "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
        "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
        "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k],
        "R" => R,
    )
    bits(x) = reinterpret(Int64, x)
    ulp(got, ref) = abs(bits(Float64(got)) - bits(Float64(ref)))

    # A single FAQ evaluation per cell must reproduce the imperative build_duo_grid
    # geometry bit-for-bit, at every subdivision level.
    for level in 0:2
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        V, F = ESD.duo_subdivide_faq(Float64, level)   # unit verts + faces, as build uses
        R = g.R
        Nc = size(F, 2)
        @test Nc == ESD.n_cells(g)
        worst = 0
        for c in 1:Nc
            i, j, k = F[1, c], F[2, c], F[3, c]
            b = bindings(V, i, j, k, R)
            worst = max(
                worst,
                ulp(ESD.eval_coeff(area, b), g.area[c]),
                ulp(ESD.eval_coeff(cart[1], b), g.cell_cart[1, c]),
                ulp(ESD.eval_coeff(cart[2], b), g.cell_cart[2, c]),
                ulp(ESD.eval_coeff(cart[3], b), g.cell_cart[3, c]),
                ulp(ESD.eval_coeff(lon, b), g.lon[c]),
                ulp(ESD.eval_coeff(lat, b), g.lat[c]),
            )
        end
        @test worst == 0   # bit-for-bit, every cell, every metric
    end

    # Physical sanity: the FAQ areas sum to the full sphere 4πR² and every cell area is
    # positive — independent confirmation the geometry is right, not just self-consistent.
    let
        g = build_duo_grid(loader = (path = "builtin://icosahedral/1", reader = "builtin_icosahedral"))
        V, F = ESD.duo_subdivide_faq(Float64, 1)
        R = g.R
        faq_area = [ESD.eval_coeff(area, bindings(V, F[1, c], F[2, c], F[3, c], R)) for c in 1:size(F, 2)]
        @test all(>(0), faq_area)
        @test sum(faq_area) ≈ 4π * R^2 rtol = 1e-12
        # lon/lat from the FAQ are the geographic image of the FAQ cartesian centroid.
        for c in 1:size(F, 2)
            b = bindings(V, F[1, c], F[2, c], F[3, c], R)
            @test ESD.eval_coeff(lat, b) ≈ asin(clamp(ESD.eval_coeff(cart[3], b) / R, -1.0, 1.0)) atol = 1e-14
            @test ESD.eval_coeff(lon, b) ≈ atan(ESD.eval_coeff(cart[2], b), ESD.eval_coeff(cart[1], b)) atol = 1e-14
        end
    end
end

@testitem "duo primal geometry FAQ: declarative doc + cross-binding byte contract" tags = [:grid, :duo, :faq, :geometry, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const ESS = EarthSciSerialization
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    FAQ_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo", "faq")

    # --- The declarative value-invention FAQ document exists, is schema-valid, and
    #     expresses the geometry chain over the cell set.
    doc_path = joinpath(FAQ_DIR, "primal_geometry.esm")
    doc = JSON.parsefile(doc_path; dicttype = Dict{String, Any})
    @test isempty(ESS.validate_schema(doc))   # validates against the bundled esm-schema.json

    model = doc["models"]["IcosahedralPrimalGeometry"]
    isets = model["index_sets"]
    for s in ("cells", "verts", "corners", "space")
        @test haskey(isets, s)
    end
    @test isets["corners"]["size"] == 3 && isets["space"]["size"] == 3

    vars = model["variables"]
    for v in ("face_vert", "vert_coord", "R", "cent_sum", "cent_unit", "cell_cart", "lon", "lat", "area")
        @test haskey(vars, v)
    end
    @test vars["face_vert"]["type"] == "parameter"
    @test vars["vert_coord"]["type"] == "parameter"
    @test vars["area"]["type"] == "state"

    eqs = model["equations"]
    lhs_name(e) = (e["lhs"] isa Dict && get(e["lhs"], "op", "") == "index") ? e["lhs"]["args"][1] : nothing
    eqfor(name) = only(filter(e -> lhs_name(e) == name, eqs))["rhs"]
    # Every output the bead names is produced by an aggregate (sum_product FAQ).
    for name in ("cell_cart", "lon", "lat", "area")
        @test eqfor(name)["op"] == "aggregate"
    end
    # The centroid sums over the three corners (the bead's "sum_product FAQ over the 3
    # vertices").
    @test haskey(eqfor("cent_sum")["ranges"], "k")            # reduction index over corners
    @test eqfor("cent_sum")["ranges"]["k"]["from"] == "corners"

    # Recursively collect the set of operator names actually used by an equation's rhs
    # (the `"op"` fields), ignoring prose (`_comment`, `description`). This checks the
    # math, not the documentation text.
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
    # lon/lat use atan2/asin of the unit centroid; the area's L'Huilier form uses the
    # great-circle arcs (acos) and the spherical-excess leaves (tan/atan/sqrt).
    @test "atan2" in ops_of(eqfor("lon"))
    @test "asin" in ops_of(eqfor("lat"))
    @test "acos" in ops_of(eqfor("arc_a"))
    let a = ops_of(eqfor("area"))
        @test "tan" in a && "atan" in a && "sqrt" in a
    end
    # Deterministic formulas, no clipping (acceptance contract): no clamp/clip/max OP in
    # any equation, and squares are products (no `^` op) so the float ops match the
    # imperative bit-for-bit. (The prose freely discusses clamp/clipping; we check ops.)
    all_ops = reduce((a, e) -> ops_of(e["rhs"], a), eqs; init = Set{String}())
    @test isempty(intersect(all_ops, Set(["clamp", "clip", "max", "^"])))

    # --- Cross-binding byte/value contract: the pinned level-0 golden is reproduced
    #     bit-for-bit by BOTH the imperative grid and the FAQ (ESS eval_coeff).
    golden = JSON.parsefile(
        joinpath(FAQ_DIR, "fixtures", "canonical", "primal_geometry_level0.json");
        dicttype = Dict{String, Any},
    )
    g = build_duo_grid(loader = (path = "builtin://icosahedral/0", reader = "builtin_icosahedral"))
    Nc = ESD.n_cells(g)
    @test golden["n_cells"] == Nc
    @test golden["level"] == 0
    bits(x) = reinterpret(Int64, Float64(x))

    # (a) imperative grid == golden, bit-for-bit.
    for c in 1:Nc
        @test bits(golden["area"][c]) == bits(g.area[c])
        @test bits(golden["cart_x"][c]) == bits(g.cell_cart[1, c])
        @test bits(golden["cart_y"][c]) == bits(g.cell_cart[2, c])
        @test bits(golden["cart_z"][c]) == bits(g.cell_cart[3, c])
        @test bits(golden["lon"][c]) == bits(g.lon[c])
        @test bits(golden["lat"][c]) == bits(g.lat[c])
    end

    # (b) the FAQ (driven through the ESS engine) == golden, bit-for-bit — the geometry
    #     is genuinely produced by the FAQ, not just by the imperative loop.
    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nrm = mk("sqrt", mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz)))
    ux = mk("/", mx, nrm); uy = mk("/", my, nrm); uz = mk("/", mz, nrm)
    area_x(d) = mk("*", "R", d)
    lon = mk("atan2", uy, ux); lat = mk("asin", uz)
    dot_bc = mk("+", mk("*", "b1", "c1"), mk("*", "b2", "c2"), mk("*", "b3", "c3"))
    dot_ca = mk("+", mk("*", "c1", "a1"), mk("*", "c2", "a2"), mk("*", "c3", "a3"))
    dot_ab = mk("+", mk("*", "a1", "b1"), mk("*", "a2", "b2"), mk("*", "a3", "b3"))
    da = mk("acos", dot_bc); db = mk("acos", dot_ca); dc = mk("acos", dot_ab)
    s = mk("*", 0.5, mk("+", da, db, dc)); half(x) = mk("/", x, 2.0)
    t = mk("*", mk("tan", half(s)), mk("tan", half(mk("-", s, da))),
        mk("tan", half(mk("-", s, db))), mk("tan", half(mk("-", s, dc))))
    area = mk("*", mk("*", 4.0, mk("atan", mk("sqrt", t))), mk("*", "R", "R"))

    V, F = ESD.duo_subdivide_faq(Float64, 0)
    R = g.R
    for c in 1:Nc
        i, j, k = F[1, c], F[2, c], F[3, c]
        b = Dict{String, Float64}(
            "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
            "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
            "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k], "R" => R,
        )
        @test bits(ESD.eval_coeff(area, b)) == bits(golden["area"][c])
        @test bits(ESD.eval_coeff(area_x(ux), b)) == bits(golden["cart_x"][c])
        @test bits(ESD.eval_coeff(area_x(uy), b)) == bits(golden["cart_y"][c])
        @test bits(ESD.eval_coeff(area_x(uz), b)) == bits(golden["cart_z"][c])
        @test bits(ESD.eval_coeff(lon, b)) == bits(golden["lon"][c])
        @test bits(ESD.eval_coeff(lat, b)) == bits(golden["lat"][c])
    end
end
