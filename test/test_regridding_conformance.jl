# Regridding cross-binding conformance — REFERENCE BINDING (Julia) consumer.
#
# esd-47z.6 (RFC pure-io-data-loaders §8, "ESD DRAIN"). The regridding/ family ships
# per-kernel Julia evaluator tests co-located with each .esm (regridding/
# {conservative_regrid_overlap_join,bspline_regrid,point_cell_average_regrid}_conformance_test.jl,
# beads esd-9cb/esd-47z.3/.4). This file is the SHARED cross-binding layer: it
# consumes the binding-neutral fixtures + golden under tests/conformance/regridding/
# <model>/ and asserts the regridding INVARIANTS are green THROUGH THE ESS ENGINE at
# the reference binding:
#
#   • bspline      — polynomial-degree reproduction on staggered locations
#   • point        — cell mean of contributing stations + configured missing_value
#   • conservative — global conservation + partition-of-unity (CONFORMANCE_SPEC §5.8)
#
# The cross-binding disposition (which invariants are binding-independent vs which
# are engine-bound — geometry clip, value-invention front-door) is documented in
# tests/conformance/regridding/README.md.

@testsnippet RegriddingConformanceSetup begin
    using Test
    import EarthSciSerialization as ESS
    # Activate EarthSciSerializationGeometryOpsExt for the conservative spherical clip.
    import GeometryOps, GeoInterface
    const JSON3 = ESS.JSON3
    using JSON

    const REGRID_HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "regridding")
    const REPO_ROOT = joinpath(@__DIR__, "..")

    _raw(esm_relpath) = JSON3.read(read(joinpath(REPO_ROOT, esm_relpath), String))
    close_abs(a, b, atol) = abs(a - b) <= atol
    close_rel(a, b, rtol) = abs(a - b) <= rtol * max(1.0, abs(b))

    # ---- B-spline driver: F_tgt[j] = Σ_k w_k(s[j])·F_src[base[j]+k] (zero-IC surfacing).
    function bs_eval(esm, model_name, ntgt, const_arrays)
        ics = Dict("F_tgt[$j]" => 0.0 for j in 1:ntgt)
        f!, u0, p, _, vmap = ESS.build_evaluator(
            _raw(esm); model_name = model_name,
            initial_conditions = ics, const_arrays = const_arrays,
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        return [du[vmap["F_tgt[$j]"]] for j in 1:ntgt]
    end

    # ---- Point driver: count[j], F_tgt[j] from the value-invention front-door.
    function pt_eval(esm, model_name, n, ca, dx, dy, missing_value)
        ics = Dict(
            ("count[$j]" => 0.0 for j in 1:n)...,
            ("F_tgt[$j]" => 0.0 for j in 1:n)...,
        )
        f!, u0, p, _, vmap = ESS.build_evaluator(
            _raw(esm); model_name = model_name,
            initial_conditions = ics, const_arrays = ca,
            parameter_overrides = Dict("dx" => dx, "dy" => dy, "missing_value" => missing_value),
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        count = [du[vmap["count[$j]"]] for j in 1:n]
        F_tgt = [du[vmap["F_tgt[$j]"]] for j in 1:n]
        return count, F_tgt
    end

    # ---- Conservative driver: build A_ij via the spherical clip, then assemble.
    ring_matrix(ring) = Float64[ring[k][d] for k in 1:length(ring), d in 1:2]

    function cons_eval(esm, model_name, src, tgt, F_src, src_lon, src_lat, tgt_lon, tgt_lat, repi, repj, dx, dy, atol)
        nsrc = length(src); ntgt = length(tgt)
        A_ij = zeros(Float64, nsrc, ntgt)
        for i in 1:nsrc, j in 1:ntgt
            ring = ESS.intersect_polygon(src[i], tgt[j], "spherical")
            size(ring, 1) < 3 && continue
            A_ij[i, j] = ESS.polygon_area(ring, "spherical")
        end
        dst_areas = vec(sum(A_ij; dims = 1))
        src_areas = vec(sum(A_ij; dims = 2))
        ca = Dict(
            "A_ij" => A_ij, "F_src" => F_src, "dst_areas" => dst_areas,
            "src_lon" => src_lon, "src_lat" => src_lat,
            "tgt_lon" => tgt_lon, "tgt_lat" => tgt_lat,
            "src_poly_rep" => src[repi], "tgt_poly_rep" => tgt[repj],
        )
        ics = Dict(
            "narrow_phase_area" => 0.0,
            ("A_j[$j]" => 0.0 for j in 1:ntgt)...,
            ("F_tgt[$j]" => 0.0 for j in 1:ntgt)...,
        )
        f!, u0, p, _, vmap = ESS.build_evaluator(
            _raw(esm); model_name = model_name,
            initial_conditions = ics, const_arrays = ca,
            parameter_overrides = Dict("dx" => dx, "dy" => dy, "atol" => atol),
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        A_j = [du[vmap["A_j[$j]"]] for j in 1:ntgt]
        F_tgt = [du[vmap["F_tgt[$j]"]] for j in 1:ntgt]
        return A_ij, A_j, F_tgt, dst_areas, src_areas
    end
end

@testitem "regridding bspline cross-binding fixture (esd-47z.6)" setup = [RegriddingConformanceSetup] tags = [:regridding, :bspline, :conformance] begin
    dir = joinpath(REGRID_HARNESS, "bspline")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; atol = Float64(fx["tolerance"]["absolute"])
    @test isfile(joinpath(REPO_ROOT, esm))
    @test length(gold["models"]) == length(fx["models"])

    for gm in gold["models"]
        name = gm["name"]
        F_tgt_exp = Float64.(gm["F_tgt"]); ntgt = length(F_tgt_exp)
        @testset "bspline $name (polynomial reproduction)" begin
            if haskey(gm, "base")   # 1-D model
                ca = Dict(
                    "F_src" => Float64.(gm["F_src"]),
                    "base" => Float64.(gm["base"]), "s" => Float64.(gm["s"]),
                )
                got = bs_eval(esm, name, ntgt, ca)
            else                    # 2-D bilinear
                fs = gm["F_src"]; nx = length(fs); ny = length(fs[1])
                F_src = Float64[Float64(fs[i][j]) for i in 1:nx, j in 1:ny]
                ca = Dict(
                    "F_src" => F_src,
                    "base_x" => Float64.(gm["base_x"]), "base_y" => Float64.(gm["base_y"]),
                    "s_x" => Float64.(gm["s_x"]), "s_y" => Float64.(gm["s_y"]),
                )
                got = bs_eval(esm, name, ntgt, ca)
            end
            for j in 1:ntgt
                @test close_abs(got[j], F_tgt_exp[j], atol)   # engine reproduces the polynomial oracle
            end
        end
    end
end

@testitem "regridding point cross-binding fixture (esd-47z.6)" setup = [RegriddingConformanceSetup] tags = [:regridding, :point, :conformance] begin
    dir = joinpath(REGRID_HARNESS, "point")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; model = fx["esm_model_name"]; atol = Float64(fx["tolerance"]["absolute"])
    dx = Float64(fx["dx"]); dy = Float64(fx["dy"])
    mv = Float64(fx["missing_value"])
    @test isfile(joinpath(REPO_ROOT, esm))

    count_exp = Float64.(gold["count"]); F_exp = Float64.(gold["F_tgt"])
    n = length(F_exp)
    ca = Dict(
        "station_lon" => Float64.(fx["stations"]["lon"]),
        "station_lat" => Float64.(fx["stations"]["lat"]),
        "cell_lon" => Float64.(fx["cells"]["lon"]),
        "cell_lat" => Float64.(fx["cells"]["lat"]),
        "station_val" => Float64.(fx["stations"]["val"]),
        "count_denom" => Float64.(gold["count_denom"]),
    )

    count, F_tgt = pt_eval(esm, model, n, ca, dx, dy, mv)
    @testset "point cell mean + missing_value" begin
        for j in 1:n
            @test close_abs(count[j], count_exp[j], atol)
            @test close_abs(F_tgt[j], F_exp[j], atol)
        end
        # the missing fill is genuinely configurable (the empty cell tracks it)
        alt = Float64(fx["alt_missing_value"])
        _, F_alt = pt_eval(esm, model, n, ca, dx, dy, alt)
        empties = findall(==(0.0), count_exp)
        for j in empties
            @test close_abs(F_alt[j], alt, atol)
        end
    end
end

@testitem "regridding conservative cross-binding fixture (esd-47z.6)" setup = [RegriddingConformanceSetup] tags = [:regridding, :conservative, :conformance] begin
    dir = joinpath(REGRID_HARNESS, "conservative")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; model = fx["esm_model_name"]
    pou_atol = Float64(fx["tolerance"]["partition_of_unity_absolute"])
    cons_rtol = Float64(fx["tolerance"]["conservation_relative"])
    area_rtol = Float64(fx["tolerance"]["area_relative"])
    @test isfile(joinpath(REPO_ROOT, esm))

    src = [ring_matrix(r) for r in fx["src_polys"]]
    tgt = [ring_matrix(r) for r in fx["tgt_polys"]]
    F_src = Float64.(fx["F_src"])
    repi = Int(fx["narrow_phase_rep"]["src"]); repj = Int(fx["narrow_phase_rep"]["tgt"])
    A_ij, A_j, F_tgt, dst_areas, src_areas = cons_eval(
        esm, model, src, tgt, F_src,
        Float64.(fx["src_rep_lon"]), Float64.(fx["src_rep_lat"]),
        Float64.(fx["tgt_rep_lon"]), Float64.(fx["tgt_rep_lat"]),
        repi, repj, Float64(fx["dx"]), Float64(fx["dy"]), Float64(fx["atol"]),
    )
    ntgt = length(A_j)

    @testset "conservative invariants through the ESS engine" begin
        # A_j reproduces the build-once dst_areas (the buffer-gated join assembly)
        @test A_j ≈ dst_areas
        # PARTITION-OF-UNITY: Σ_i W_ij = 1, exact by construction
        for j in 1:ntgt
            w_sum = sum(A_ij[i, j] for i in 1:length(src)) / A_j[j]
            @test close_abs(w_sum, 1.0, pou_atol)
        end
        # CONSERVATION: Σ_j A_j·F_tgt = Σ_i A_i·F_src
        @test close_rel(sum(A_j .* F_tgt), sum(src_areas .* F_src), cons_rtol)
        # regression against the reference-binding golden (area tolerance, §5.8)
        F_gold = Float64.(gold["F_tgt"]); A_gold = Float64.(gold["A_j"])
        for j in 1:ntgt
            @test close_rel(F_tgt[j], F_gold[j], area_rtol)
            @test close_rel(A_j[j], A_gold[j], area_rtol)
        end
        @test close_rel(sum(A_j .* F_tgt), Float64(gold["target_mass"]), area_rtol)
    end
end
