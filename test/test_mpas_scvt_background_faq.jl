# MPAS-SCVT background quadrature mesh + density sampling as a declarative CONST
# integrand FAQ (esd-e5m.1 / D1).
#
# This is the FIXED point set the SCVT Lloyd step (esd-e5m.2 / D2) integrates
# over. One iteration forms, per generator g, the density-weighted centroid
#   c_g = (sum_{p->g} rho_p dA_p x_p) / (sum_{p->g} rho_p dA_p)
# over the background quadrature points p assigned to g, then moves g to c_g
# (re-normalized). D1 declares the static, generator-independent half: the
# quadrature point set + the per-point density-weighted integrands the step sums.
#
# The point set is DUO-subdivision-based, reusing the esd-heg.3 subdivision idiom
# (subdivide_fold.esm / esd-heg.4 over the base seed): each quadrature point IS a
# DUO icosahedral primal cell, its position bg_coord the unit-sphere cell centroid
# and its weight `area` the spherical-triangle dA (primal_geometry.esm / esd-heg.6).
# The FAQ's new declarative content is the integrand: bg_mass = rho * area (the
# centroid denominator) and bg_moment = bg_mass * bg_coord (the numerator).
# Density sampling: rho_p = f(bg_coord_p) through the single-pathway evaluator
# (eval_coeff), the same mechanism as the expression-IC rule; rho == 1 is the
# uniform-density regression.
#
# Per AGENTS.md (single pathway), every value is produced through the landed ESS
# engine via ESD.eval_coeff — no shadow evaluator. The scalar per-cell product
# ASTs driven here are the scalar instances of background_quadrature.esm's indexed
# sum_product aggregates (q-1 / q-2), whose indexed structure is checked
# separately against the schema-valid document.

@testitem "scvt background quadrature: integrand via FAQ + physical invariants" tags = [:grid, :mpas, :scvt, :faq, :quadrature] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    using LinearAlgebra: norm

    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

    # Centroid AST over the three corner unit vectors (= primal_geometry ux/uy/uz):
    # u = (a+b+c) / |a+b+c|, squares as products. The unit-sphere quadrature point.
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nrm = mk("sqrt", mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz)))
    ucomp = (mk("/", mx, nrm), mk("/", my, nrm), mk("/", mz, nrm))

    # Density sampling AST: rho(x,y,z) = 2 + z (smooth, latitude-graded, positive).
    rho_ast = mk("+", 2.0, "z")
    # FAQ integrand product ASTs (the scalar instances of q-1 / q-2).
    mass_ast = mk("*", "rho", "area")        # bg_mass   = rho * area
    moment_ast = mk("*", "mass", "coord")    # bg_moment = bg_mass * bg_coord

    cb(V, i, j, k) = Dict{String, Float64}(
        "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
        "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
        "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k],
    )
    bits(x) = reinterpret(Int64, Float64(x))

    for level in 0:2
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        V, F = ESD.duo_subdivide_faq(Float64, level)
        R = g.R
        Nc = ESD.n_cells(g)

        coord = Matrix{Float64}(undef, 3, Nc)
        area = Vector{Float64}(undef, Nc)
        for c in 1:Nc
            b = cb(V, F[1, c], F[2, c], F[3, c])
            coord[1, c] = ESD.eval_coeff(ucomp[1], b)
            coord[2, c] = ESD.eval_coeff(ucomp[2], b)
            coord[3, c] = ESD.eval_coeff(ucomp[3], b)
            area[c] = g.area[c]
        end

        # Density sampling through the single pathway.
        rho1 = ones(Float64, Nc)
        rhos = [ESD.eval_coeff(rho_ast, Dict{String, Float64}("z" => coord[3, c])) for c in 1:Nc]

        # FAQ integrands via eval_coeff, for each density.
        function faq_integrands(rho)
            mass = [ESD.eval_coeff(mass_ast, Dict{String, Float64}("rho" => rho[c], "area" => area[c])) for c in 1:Nc]
            mom = Matrix{Float64}(undef, 3, Nc)
            for c in 1:Nc, d in 1:3
                mom[d, c] = ESD.eval_coeff(moment_ast, Dict{String, Float64}("mass" => mass[c], "coord" => coord[d, c]))
            end
            return mass, mom
        end
        mass1, mom1 = faq_integrands(rho1)
        masss, moms = faq_integrands(rhos)

        # (A) rho == 1 reduces the integrand to the bare quadrature weight, bit-for-bit:
        #     bg_mass == area and bg_moment == area * centroid.
        for c in 1:Nc
            @test bits(mass1[c]) == bits(area[c])
            for d in 1:3
                @test bits(mom1[d, c]) == bits(area[c] * coord[d, c])
            end
        end

        # (B) Density sampling is exactly rho = 2 + z, and always positive.
        for c in 1:Nc
            @test bits(rhos[c]) == bits(2.0 + coord[3, c])
        end
        @test all(>(0), rhos)

        # (C) The quadrature weights tile the sphere: sum of the rho==1 measure is
        #     the full sphere area 4 pi R^2 (and every weight is positive).
        @test all(>(0), mass1)
        @test sum(mass1) ≈ 4π * R^2 rtol = 1e-12

        # (D) The density-weighted centroid (the quantity the step divides to get c_g):
        #     uniform density -> sphere centre (origin); the rho = 2 + z density ->
        #     shifts toward the dense (+z) hemisphere, to the analytic value 1/6
        #     ( = integral (2+z) z dA / integral (2+z) dA over the unit sphere ).
        cen(mass, mom) = [sum(@view mom[d, :]) / sum(mass) for d in 1:3]
        c_uniform = cen(mass1, mom1)
        c_sampled = cen(masss, moms)
        @test norm(c_uniform) < 1e-9   # ~0 (origin) by icosahedral symmetry
        @test c_sampled[3] > 0.1             # shifted north toward the dense hemisphere
        @test c_sampled[3] ≈ 1 / 6 atol = 2e-2
        @test abs(c_sampled[1]) < 1e-9 && abs(c_sampled[2]) < 1e-9
    end
end

@testitem "scvt background quadrature FAQ: declarative doc + canonical byte contract" tags = [:grid, :mpas, :scvt, :faq, :quadrature, :conformance] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using EarthSciSerialization
    const ESS = EarthSciSerialization
    using JSON

    REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    SCVT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt")

    # --- The declarative document exists, is schema-valid, and expresses the
    #     density-weighted integrand the SCVT step sums.
    doc_path = joinpath(SCVT_DIR, "background_quadrature.esm")
    doc = JSON.parsefile(doc_path; dicttype = Dict{String, Any})
    @test isempty(ESS.validate_schema(doc))   # validates against the bundled esm-schema.json

    model = doc["models"]["ScvtBackgroundQuadrature"]
    isets = model["index_sets"]
    @test isets["cells"]["kind"] == "interval"
    @test isets["space"]["size"] == 3

    vars = model["variables"]
    # The quadrature point set + sampled density are CONST inputs (parameters);
    # the integrands the step sums are derived state.
    @test vars["bg_coord"]["type"] == "parameter" && vars["bg_coord"]["shape"] == ["cells", "space"]
    @test vars["area"]["type"] == "parameter" && vars["area"]["shape"] == ["cells"]
    @test vars["rho"]["type"] == "parameter" && vars["rho"]["shape"] == ["cells"]
    @test vars["bg_mass"]["type"] == "state" && vars["bg_mass"]["shape"] == ["cells"]
    @test vars["bg_moment"]["type"] == "state" && vars["bg_moment"]["shape"] == ["cells", "space"]

    eqs = model["equations"]
    lhs_name(e) = (e["lhs"] isa Dict && get(e["lhs"], "op", "") == "index") ? e["lhs"]["args"][1] : nothing
    eqfor(name) = only(filter(e -> lhs_name(e) == name, eqs))["rhs"]
    # bg_mass = rho * area: an elementwise sum_product aggregate over cells.
    massrhs = eqfor("bg_mass")
    @test massrhs["op"] == "aggregate"
    @test get(massrhs, "semiring", "sum_product") == "sum_product"
    @test massrhs["output_idx"] == ["c"]
    @test massrhs["ranges"]["c"]["from"] == "cells"
    @test massrhs["expr"]["op"] == "*"
    # bg_moment = bg_mass * bg_coord: aggregate over cell x space.
    momrhs = eqfor("bg_moment")
    @test momrhs["op"] == "aggregate"
    @test momrhs["output_idx"] == ["c", "d"]
    @test momrhs["ranges"]["d"]["from"] == "space"
    @test momrhs["expr"]["op"] == "*"

    # Deterministic products only — no clipping, no `^` (squares-as-products bit
    # stability); the only ops are the product and the index reads.
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
    @test isempty(intersect(all_ops, Set(["clamp", "clip", "max", "min", "^"])))
    @test all_ops == Set(["aggregate", "*", "index"])

    # --- Cross-binding byte/value contract: the pinned goldens are reproduced
    #     bit-for-bit by the FAQ (ESS eval_coeff) AND the imperative grid (area).
    mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
    mx = mk("+", "a1", "b1", "c1"); my = mk("+", "a2", "b2", "c2"); mz = mk("+", "a3", "b3", "c3")
    nrm = mk("sqrt", mk("+", mk("*", mx, mx), mk("*", my, my), mk("*", mz, mz)))
    ucomp = (mk("/", mx, nrm), mk("/", my, nrm), mk("/", mz, nrm))
    rho_ast = mk("+", 2.0, "z")
    mass_ast = mk("*", "rho", "area")
    moment_ast = mk("*", "mass", "coord")
    cb(V, i, j, k) = Dict{String, Float64}(
        "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
        "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
        "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k],
    )
    bits(x) = reinterpret(Int64, Float64(x))

    for level in (0, 1)
        golden = JSON.parsefile(
            joinpath(SCVT_DIR, "fixtures", "canonical", "background_quadrature_level$(level).json");
            dicttype = Dict{String, Any},
        )
        g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
        V, F = ESD.duo_subdivide_faq(Float64, level)
        Nc = ESD.n_cells(g)
        @test golden["n_cells"] == Nc
        @test golden["level"] == level

        for c in 1:Nc
            b = cb(V, F[1, c], F[2, c], F[3, c])
            cx = ESD.eval_coeff(ucomp[1], b); cy = ESD.eval_coeff(ucomp[2], b); cz = ESD.eval_coeff(ucomp[3], b)
            area = g.area[c]

            # quadrature point + weight reproduce the golden bit-for-bit (FAQ + imperative).
            @test bits(cx) == bits(golden["bg_coord_x"][c])
            @test bits(cy) == bits(golden["bg_coord_y"][c])
            @test bits(cz) == bits(golden["bg_coord_z"][c])
            @test bits(area) == bits(golden["area"][c])

            # uniform-density integrand (rho == 1).
            m1 = ESD.eval_coeff(mass_ast, Dict{String, Float64}("rho" => 1.0, "area" => area))
            @test bits(m1) == bits(golden["bg_mass_uniform"][c])
            for (d, key) in zip(1:3, ("bg_moment_uniform_x", "bg_moment_uniform_y", "bg_moment_uniform_z"))
                cd = (cx, cy, cz)[d]
                @test bits(ESD.eval_coeff(moment_ast, Dict{String, Float64}("mass" => m1, "coord" => cd))) == bits(golden[key][c])
            end

            # sampled-density integrand (rho = 2 + z).
            rs = ESD.eval_coeff(rho_ast, Dict{String, Float64}("z" => cz))
            @test bits(rs) == bits(golden["rho_sampled"][c])
            ms = ESD.eval_coeff(mass_ast, Dict{String, Float64}("rho" => rs, "area" => area))
            @test bits(ms) == bits(golden["bg_mass_sampled"][c])
            for (d, key) in zip(1:3, ("bg_moment_sampled_x", "bg_moment_sampled_y", "bg_moment_sampled_z"))
                cd = (cx, cy, cz)[d]
                @test bits(ESD.eval_coeff(moment_ast, Dict{String, Float64}("mass" => ms, "coord" => cd))) == bits(golden[key][c])
            end
        end
    end
end
