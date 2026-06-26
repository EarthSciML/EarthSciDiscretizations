# Reduction cross-binding conformance — REFERENCE BINDING (Julia) consumer.
#
# Bead esd-khv (campfire-e2e C4 design spike 2026-06-26). The reduction/ family
# ships the per-kernel Julia evaluator test co-located with the .esm
# (reduction/lev_min_surface_reduce_conformance_test.jl) with the rich end-to-end
# assertions; this file is the SHARED cross-binding layer: it consumes the
# binding-neutral fixtures + golden under tests/conformance/reduction/<model>/ and
# asserts the reduction reproduces the analytic oracle THROUGH THE ESS ENGINE at
# the reference binding.
#
#   • lev_min — F_surf[x,y] = F_3d[x,y, argmin_lev lev_coord]; the ground_surface
#               interface's dimension_mapping.constraints.lev = {value:min}.
#
# Like bspline (regridding), lev_min is the family's BINDING-INDEPENDENT,
# closed-form tier: no value-invention front-door, no geometry clip. The
# cross-binding disposition is documented in tests/conformance/reduction/README.md.

@testsnippet ReductionConformanceSetup begin
    using Test
    import EarthSciSerialization as ESS
    const JSON3 = ESS.JSON3
    using JSON

    const REDUCTION_HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "reduction")
    const REPO_ROOT = joinpath(@__DIR__, "..")

    _raw(esm_relpath) = JSON3.read(read(joinpath(REPO_ROOT, esm_relpath), String))
    close_abs(a, b, atol) = abs(a - b) <= atol

    # Reconstruct the (nx,ny,nlev) source field from the nested fixture list
    # F_3d[x][y][lev], and drive the reduction through the ESS engine. du = f!(u0)
    # at the zero IC IS the surfaced field (constant-RHS D-equation).
    function lev_min_eval(esm, model_name, F_3d_nested, lev_coord)
        nx = length(F_3d_nested); ny = length(F_3d_nested[1]); nlev = length(F_3d_nested[1][1])
        F_3d = Array{Float64}(undef, nx, ny, nlev)
        for a in 1:nx, b in 1:ny, l in 1:nlev
            F_3d[a, b, l] = Float64(F_3d_nested[a][b][l])
        end
        ca = Dict("F_3d" => F_3d, "lev_coord" => Float64.(lev_coord))
        ics = Dict("F_surf[$a,$b]" => 0.0 for a in 1:nx, b in 1:ny)
        f!, u0, p, _, vmap = ESS.build_evaluator(
            _raw(esm); model_name = model_name,
            initial_conditions = ics, const_arrays = ca,
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        return [du[vmap["F_surf[$a,$b]"]] for a in 1:nx, b in 1:ny], F_3d
    end
end

@testitem "reduction lev_min cross-binding fixture (esd-khv)" setup = [ReductionConformanceSetup] tags = [:reduction, :lev_min, :conformance] begin
    dir = joinpath(REDUCTION_HARNESS, "lev_min")
    fx = JSON.parsefile(joinpath(dir, "fixtures.json"))
    gold = JSON.parsefile(joinpath(dir, "golden.json"))
    esm = fx["esm"]; model = fx["esm_model_name"]; atol = Float64(fx["tolerance"]["absolute"])
    @test isfile(joinpath(REPO_ROOT, esm))

    F_surf, F_3d = lev_min_eval(esm, model, fx["F_3d"], fx["lev_coord"])
    nx = size(F_surf, 1); ny = size(F_surf, 2)
    F_exp = gold["F_surf"]                       # nested [x][y]
    lmin = Int(gold["argmin_level"])

    @testset "lev=min surface == analytic argmin slice (through the ESS engine)" begin
        for a in 1:nx, b in 1:ny
            @test close_abs(F_surf[a, b], Float64(F_exp[a][b]), atol)   # engine reproduces the oracle
        end
        # the golden IS the F_3d slice at the minimum-coordinate level
        for a in 1:nx, b in 1:ny
            @test F_surf[a, b] == F_3d[a, b, lmin]
        end
        # the argmin search is LOAD-BEARING: the canonical fixture's minimum is at
        # an interior level, so neither end-slice equals the surfaced field.
        nlev = size(F_3d, 3)
        @test 1 < lmin < nlev                                  # interior minimum
        @test [F_surf[a, b] for a in 1:nx, b in 1:ny] != F_3d[:, :, 1]
        @test [F_surf[a, b] for a in 1:nx, b in 1:ny] != F_3d[:, :, nlev]
    end
end
