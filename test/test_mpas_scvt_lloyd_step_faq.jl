# MPAS-SCVT one-iteration Lloyd STEP as a declarative value-invention FAQ
# (esd-e5m.2 / D2), evaluated end-to-end by the EarthSciSerialization
# value-invention front-door (materialize_value_invention).
#
# ONE Lloyd / SCVT iteration over the FIXED background quadrature point set
# declared by background_quadrature.esm (D1): per generator g it forms the
# density-weighted centroid c_g = (sum_{p->g} rho_p dA_p x_p)/(sum_{p->g} rho_p dA_p)
# over the background points p assigned to g, then moves g to c_g. The three
# STEP pieces are all declarative through the front-door:
#   (1) ASSIGN   assign[c] = argmin_g dist2(bg_coord[c], gen[g])  (E1, ess-os1)
#   (2) CENTROID den[g] = sum_{c->g} bg_mass[c];
#               num_{x,y,z}[g] = sum_{c->g} bg_mass[c]*bg_coord[c,{1,2,3}] (E2, ess-2u5)
#   (3) MOVE     centroid_{x,y,z}[g] = num_{x,y,z}[g] / den[g]
# A grouped semiring reduction emits a single output index, so the 3-D moment is
# three scalar reductions (num_x/num_y/num_z), not one [generators,space] aggregate.
# The final sphere RE-PROJECTION gen_next = R*centroid/|centroid| needs a sqrt and
# is the host loop's RHS-only step (LLOYD_STEP_CONTRACT.md), not a value-invention
# buffer — the build-time relational front-door carries no transcendental geometry.
#
# Per AGENTS.md (single pathway), every value is produced through the landed ESS
# engine: the quadrature point set + density-weighted measure come from the D1 FAQ
# (eval_coeff of the same ASTs as background_quadrature.esm), and the SCVT step is
# the end-to-end output of materialize_value_invention driving lloyd_step.esm.

@testitem "scvt lloyd step: assign->centroid->move via the value-invention front-door" tags = [:grid, :mpas, :scvt, :faq, :lloyd] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const ESS = EarthSciSerialization
    using JSON
    using LinearAlgebra: norm

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    SCVT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt")
    ESM_PATH = joinpath(SCVT_DIR, "lloyd_step.esm")

    load_model() = ESS._select_model_json(
        JSON.parsefile(ESM_PATH; dicttype = Dict{String, Any}), "ScvtLloydStep"
    )
    run_step(mj, bg_coord, bg_mass, gen) = ESS.materialize_value_invention(
        mj, Dict("bg_coord" => bg_coord, "bg_mass" => bg_mass, "gen" => gen),
        Dict{String, Float64}()
    )

    # ── (A) Planar regression: the 3-D step reduces EXACTLY to the landed ESS E2
    #        fixture (nearest_generator_centroid.esm). Generators on the x-axis at
    #        0,1,2; points at 0,0.75,1.25,2.0 (exact dyadics, no ties), weights
    #        bg_mass = [1,1,3,4]. Generator 2 owns points 2,3 → centroid_x 4.5/4 =
    #        1.125 (moved from 1.0); the y,z components are identically 0.
    @testset "planar regression: reduces to the landed E2 centroid (bit-exact)" begin
        mj = load_model()
        mj["index_sets"]["cells"]["size"] = 4
        mj["index_sets"]["generators"]["size"] = 3
        bg_coord = Float64[0 0 0; 0.75 0 0; 1.25 0 0; 2.0 0 0]
        bg_mass = Float64[1, 1, 3, 4]
        gen = Float64[0 0 0; 1 0 0; 2 0 0]
        vi = run_step(mj, bg_coord, bg_mass, gen)
        @test vi.assignments["assign"] == [1, 2, 2, 3]
        @test vi.groups["den"] == [1.0, 4.0, 4.0]
        @test vi.groups["num_x"] == [0.0, 4.5, 8.0]
        @test vi.groups["centroid_x"] == [0.0, 1.125, 2.0]   # == landed E2 result
        @test vi.groups["centroid_y"] == [0.0, 0.0, 0.0]
        @test vi.groups["centroid_z"] == [0.0, 0.0, 0.0]
        # All step buffers are value-invention (materialised at setup, off the ODE).
        @test vi.vi_var_names == Set(
            [
                "assign", "den", "num_x", "num_y", "num_z",
                "centroid_x", "centroid_y", "centroid_z",
            ]
        )
        # Pure function of inputs — a fresh build is identical.
        mj2 = load_model()
        mj2["index_sets"]["cells"]["size"] = 4
        mj2["index_sets"]["generators"]["size"] = 3
        @test run_step(mj2, bg_coord, bg_mass, gen).groups["centroid_x"] == [0.0, 1.125, 2.0]
    end

    # Geometry shared by the spherical cases: the unit-sphere DUO primal-cell
    # centroid (= background_quadrature.esm bg_coord) and the icosahedral-vertex
    # generators, both through eval_coeff (the single pathway).
    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nrm = mk("sqrt", mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz)))
    ucomp = (mk("/", mx, nrm), mk("/", my, nrm), mk("/", mz, nrm))
    rho_ast = mk("+", 2.0, "z")
    mass_ast = mk("*", "rho", "area")
    cb(V, i, j, k) = Dict{String, Float64}(
        "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
        "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
        "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k]
    )

    function level_geometry(level)
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        V, F = ESD.duo_subdivide_faq(Float64, level)
        Nc = ESD.n_cells(g)
        bg = Matrix{Float64}(undef, Nc, 3)
        area = Vector{Float64}(undef, Nc)
        for c in 1:Nc
            b = cb(V, F[1, c], F[2, c], F[3, c])
            for d in 1:3
                bg[c, d] = ESD.eval_coeff(ucomp[d], b)
            end
            area[c] = g.area[c]
        end
        return g, bg, area
    end
    # The 12 unit icosahedral vertices (level-0) are the generators in every case.
    V0, _ = ESD.duo_subdivide_faq(Float64, 0)
    GEN = permutedims(V0)   # (12 x 3)

    # Independent nearest-generator argmin with the §5.7 smallest-id tie-break
    # (ascending g + strict `<` keeps the first/smallest at a tie), evaluating the
    # squared distance with the SAME left-associated arithmetic as lloyd_step.esm.
    function indep_assign(bg, gen)
        Nc = size(bg, 1); G = size(gen, 1); a = zeros(Int, Nc)
        for c in 1:Nc
            best = Inf; barg = 0
            for gg in 1:G
                dx = bg[c, 1] - gen[gg, 1]; dy = bg[c, 2] - gen[gg, 2]; dz = bg[c, 3] - gen[gg, 3]
                d2 = (dx * dx) + (dy * dy) + (dz * dz)
                d2 < best && (best = d2; barg = gg)
            end
            a[c] = barg
        end
        return a
    end

    # ── (B) Level-0 spherical composition: D1's 20 background quadrature points x
    #        12 icosahedral generators, driven end-to-end through the front-door.
    @testset "level-0 composition: assignment, mass conservation, centroid" begin
        mj = load_model()
        _, bg, area = level_geometry(0)
        @test size(bg, 1) == 20 && size(GEN, 1) == 12

        for (rho, label) in (
                (ones(20), "uniform"),
                ([ESD.eval_coeff(rho_ast, Dict("z" => bg[c, 3])) for c in 1:20], "sampled"),
            )
            bg_mass = [ESD.eval_coeff(mass_ast, Dict("rho" => rho[c], "area" => area[c])) for c in 1:20]
            vi = run_step(mj, bg, bg_mass, GEN)
            assign = vi.assignments["assign"]

            # E1: the front-door assignment is exactly the independent argmin (the
            # §5.7 smallest-id tie-break) — a valid partition into 1:12.
            @test assign == indep_assign(bg, GEN)
            @test all(1 .<= assign .<= 12)

            # E2: mass conservation — the grouped reduction partitions ALL the mass.
            @test sum(vi.groups["den"]) ≈ sum(bg_mass) rtol = 1.0e-12
            # Independent recompute of den / num / centroid from assign + factors.
            for g in 1:12
                cells_g = findall(==(g), assign)
                den_g = isempty(cells_g) ? 0.0 : sum(bg_mass[c] for c in cells_g)
                @test vi.groups["den"][g] ≈ den_g rtol = 1.0e-12
                for (axis, key) in ((1, "num_x"), (2, "num_y"), (3, "num_z"))
                    num_g = isempty(cells_g) ? 0.0 : sum(bg_mass[c] * bg[c, axis] for c in cells_g)
                    @test vi.groups[key][g] ≈ num_g atol = 1.0e-3 * (sum(bg_mass) + 1)
                end
                cen = (vi.groups["centroid_x"][g], vi.groups["centroid_y"][g], vi.groups["centroid_z"][g])
                if den_g == 0.0
                    # Unattended generator: empty-group 0̄ denominator → NaN centroid.
                    @test all(isnan, cen)
                else
                    @test vi.groups["centroid_x"][g] ≈ vi.groups["num_x"][g] / vi.groups["den"][g] rtol = 1.0e-12
                    # The density-weighted centroid is a convex combination of
                    # unit-sphere points → inside the unit ball.
                    @test norm(collect(cen)) <= 1 + 1.0e-12
                end
            end
        end

        # Density drives the result: rho = 2 + z reweights the per-generator mass
        # away from the uniform measure (the variable-resolution mechanism).
        m_u = area
        m_s = [ESD.eval_coeff(mass_ast, Dict("rho" => 2.0 + bg[c, 3], "area" => area[c])) for c in 1:20]
        @test run_step(mj, bg, m_u, GEN).groups["den"] != run_step(mj, bg, m_s, GEN).groups["den"]
    end

    # ── (C) CVT fixed-point regression: the icosahedral generators are a CENTROIDAL
    #        Voronoi fixed point, so under uniform density the one-step move shrinks
    #        as the background quadrature refines, and a background finer than the
    #        generator set attends every generator (no empty groups).
    @testset "CVT regression: refinement attends every generator + shrinks the move" begin
        mj = load_model()
        disp = Float64[]
        for level in (1, 2)
            _, bg, area = level_geometry(level)
            mj_l = load_model()
            mj_l["index_sets"]["cells"]["size"] = size(bg, 1)
            vi = run_step(mj_l, bg, copy(area), GEN)
            # Every generator attended once the background is finer than the seeds.
            @test all(!=(0.0), vi.groups["den"])
            @test !any(isnan, vi.groups["centroid_x"])
            # Max displacement of the projected centroid from its own generator.
            md = 0.0
            for g in 1:12
                c = [vi.groups["centroid_x"][g], vi.groups["centroid_y"][g], vi.groups["centroid_z"][g]]
                md = max(md, norm(c ./ norm(c) .- GEN[g, :]))
            end
            push!(disp, md)
        end
        @test disp[1] < 0.2          # near a fixed point already at level 1
        @test disp[2] < disp[1]      # refinement shrinks the Lloyd step (CVT convergence)
    end
end

@testitem "scvt lloyd step FAQ: declarative doc + canonical byte contract" tags = [:grid, :mpas, :scvt, :faq, :lloyd, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const ESS = EarthSciSerialization
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    SCVT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt")

    # --- The declarative document exists, is schema-valid, and expresses the
    #     assign -> grouped centroid -> elementwise move SCVT step.
    doc_path = joinpath(SCVT_DIR, "lloyd_step.esm")
    doc = JSON.parsefile(doc_path; dicttype = Dict{String, Any})
    @test isempty(ESS.validate_schema(doc))

    model = doc["models"]["ScvtLloydStep"]
    isets = model["index_sets"]
    @test isets["cells"]["size"] == 20
    @test isets["generators"]["size"] == 12
    @test isets["space"]["size"] == 3

    vars = model["variables"]
    @test vars["bg_coord"]["type"] == "parameter" && vars["bg_coord"]["shape"] == ["cells", "space"]
    @test vars["bg_mass"]["type"] == "parameter" && vars["bg_mass"]["shape"] == ["cells"]
    @test vars["gen"]["type"] == "parameter" && vars["gen"]["shape"] == ["generators", "space"]
    @test vars["assign"]["type"] == "state" && vars["assign"]["shape"] == ["cells"]
    for v in ("den", "num_x", "num_y", "num_z", "centroid_x", "centroid_y", "centroid_z")
        @test vars[v]["type"] == "state" && vars[v]["shape"] == ["generators"]
    end

    eqs = model["equations"]
    lhs_name(e) = (e["lhs"] isa Dict && get(e["lhs"], "op", "") == "index") ? e["lhs"]["args"][1] : nothing
    eqfor(name) = only(filter(e -> lhs_name(e) == name, eqs))["rhs"]
    # (1) ASSIGN: an argmin arg-witness over the generators candidate range.
    arhs = eqfor("assign")
    @test arhs["op"] == "aggregate" && arhs["output_idx"] == ["c"]
    @test arhs["expr"]["op"] == "argmin"
    @test arhs["expr"]["arg"] == "g"
    @test arhs["expr"]["ranges"]["g"]["from"] == "generators"
    # (2) CENTROID: grouped sum_product reductions keyed on `assign`.
    for name in ("den", "num_x", "num_y", "num_z")
        rhs = eqfor(name)
        @test rhs["op"] == "aggregate" && rhs["output_idx"] == ["g"]
        @test get(rhs, "semiring", "sum_product") == "sum_product"
        @test rhs["ranges"]["g"]["from"] == "generators" && rhs["ranges"]["c"]["from"] == "cells"
        @test rhs["join"] == [Dict("on" => [["assign", "g"]])]
    end
    # (3) MOVE: an elementwise `/` derived buffer, no join, no contraction.
    for name in ("centroid_x", "centroid_y", "centroid_z")
        rhs = eqfor(name)
        @test rhs["op"] == "aggregate" && rhs["output_idx"] == ["g"]
        @test rhs["expr"]["op"] == "/"
        @test !haskey(rhs, "join")
        @test collect(keys(rhs["ranges"])) == ["g"]
    end
    # Deterministic value-invention op set only (relational + arithmetic; no sqrt /
    # trig — the sphere projection is the host loop's step, not a VI buffer).
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
    all_ops = reduce((a, e) -> ops_of(e["rhs"], a), eqs; init = Set{String}())
    @test all_ops == Set(["aggregate", "argmin", "+", "-", "*", "/", "index"])
    @test isempty(intersect(all_ops, Set(["sqrt", "^", "exp", "sin", "cos", "clamp"])))

    # --- Cross-binding contract: the canonical level-0 golden is reproduced by the
    #     front-door — `assign` byte-identical (§5.7), the grouped/derived floats
    #     bit-for-bit by the Julia reference (tolerance contract across bindings).
    golden = JSON.parsefile(
        joinpath(SCVT_DIR, "fixtures", "canonical", "lloyd_step_level0.json");
        dicttype = Dict{String, Any}
    )
    bits(x) = reinterpret(Int64, Float64(x))

    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nrm = mk("sqrt", mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz)))
    ucomp = (mk("/", mx, nrm), mk("/", my, nrm), mk("/", mz, nrm))
    rho_ast = mk("+", 2.0, "z"); mass_ast = mk("*", "rho", "area")
    cb(V, i, j, k) = Dict{String, Float64}(
        "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
        "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
        "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k]
    )

    g = build_duo_grid(loader = (path = "builtin://icosahedral/0", reader = "builtin_icosahedral"))
    V, F = ESD.duo_subdivide_faq(Float64, 0)
    Nc = ESD.n_cells(g)
    @test golden["n_cells"] == Nc && golden["n_generators"] == 12
    @test bits(golden["R"]) == bits(g.R)

    bg = Matrix{Float64}(undef, Nc, 3); area = Vector{Float64}(undef, Nc)
    for c in 1:Nc
        b = cb(V, F[1, c], F[2, c], F[3, c])
        for d in 1:3
            bg[c, d] = ESD.eval_coeff(ucomp[d], b)
        end
        area[c] = g.area[c]
        # bg_coord / area reproduce the golden bit-for-bit (the D1 quadrature).
        @test bits(bg[c, 1]) == bits(golden["bg_coord_x"][c])
        @test bits(bg[c, 2]) == bits(golden["bg_coord_y"][c])
        @test bits(bg[c, 3]) == bits(golden["bg_coord_z"][c])
        @test bits(area[c]) == bits(golden["area"][c])
    end
    gen = permutedims(V)
    for v in 1:12
        @test bits(gen[v, 1]) == bits(golden["gen_x"][v])
        @test bits(gen[v, 2]) == bits(golden["gen_y"][v])
        @test bits(gen[v, 3]) == bits(golden["gen_z"][v])
    end

    mj = ESS._select_model_json(doc, "ScvtLloydStep")
    # A golden value of `null` marks an unattended generator (NaN centroid).
    eq_or_nan(actual, want) = want === nothing ? isnan(actual) : bits(actual) == bits(Float64(want))

    for (key, rho) in (
            ("uniform", ones(Nc)),
            ("sampled", [ESD.eval_coeff(rho_ast, Dict("z" => bg[c, 3])) for c in 1:Nc]),
        )
        bg_mass = [ESD.eval_coeff(mass_ast, Dict("rho" => rho[c], "area" => area[c])) for c in 1:Nc]
        vi = ESS.materialize_value_invention(
            mj,
            Dict("bg_coord" => bg, "bg_mass" => bg_mass, "gen" => gen), Dict{String, Float64}()
        )
        block = golden[key]
        # assign — byte-identical integer buffer (§5.7), shared by both densities.
        @test vi.assignments["assign"] == Int.(golden["assign"])
        for c in 1:Nc
            @test bits(bg_mass[c]) == bits(block["bg_mass"][c])
        end
        for v in 1:12
            @test bits(vi.groups["den"][v]) == bits(block["den"][v])
            for nk in ("num_x", "num_y", "num_z")
                @test bits(vi.groups[nk][v]) == bits(block[nk][v])
            end
            for ck in ("centroid_x", "centroid_y", "centroid_z")
                @test eq_or_nan(vi.groups[ck][v], block[ck][v])
            end
        end
    end
end
