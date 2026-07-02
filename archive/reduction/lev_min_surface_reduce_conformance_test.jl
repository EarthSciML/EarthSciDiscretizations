# lev=min 3D->surface reduction — Julia evaluator conformance.
#
# ESD owns this declarative reduction PROGRAM (lev_min_surface_reduce.esm,
# co-located here); the EarthSciSerialization engine (an ESD dependency)
# EVALUATES it. The ground_surface interface's
# `dimension_mapping.constraints.lev = {value:min}` (campfire-e2e C4 design spike
# 2026-06-26, bead esd-khv): collapse a 3-D field F_3d[x,y,lev] to the 2-D
# surface field F_surf[x,y] at the level where the vertical coordinate is MINIMUM.
#
# This is the third tier of the reprojection/regridding rule family's
# cross-binding spectrum: a binding-INDEPENDENT, closed-form (analytic) reduction
# like bspline_regrid — NO value-invention front-door (unlike point), NO geometry
# clip (unlike conservative). The reduction is the value-at-argmin idiom expressed
# PURELY over existing aggregate vocab — an outer `sum_product` that carries the
# spatial output and contracts lev, whose body keeps a level iff its coordinate
# equals the INLINE `min_sum` reduction MIN_k lev_coord[k]:
#
#   F_surf[x,y] = SUM_lev ifelse(lev_coord[lev] == MIN_k lev_coord[k], F_3d[x,y,lev], 0)
#
# NO new engine op — the declarative-or-fail gate (lev=min reduces to `aggregate`)
# PASSES. `min` returns a lev_coord element bit-for-bit, so the `==` selects
# exactly one level byte-exactly and the surfaced value is a verbatim F_3d gather:
# the reduction is 0-ULP across bindings. The single state F_surf is a zero-IC
# constant-RHS D-equation, so du = f!(u0) IS the surfaced field (the
# regridding-assembly precedent); F_3d / lev_coord are host const_arrays.

@testsnippet LevMinReduceSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach the
    # raw document through the dependency so it is the exact type build_evaluator
    # expects.
    const JSON3 = ESS.JSON3

    const _LEVMIN_FIXTURE = joinpath(@__DIR__, "lev_min_surface_reduce.esm")

    # A 2x2x3 source field whose last decimal digit IS the level index, so the
    # surfaced slice is human-readable: F_3d[a,b,lev] = 100a + 10b + lev.
    _make_F3d() = Float64[100a + 10b + l for a in 1:2, b in 1:2, l in 1:3]

    # Drive the reduction through build_evaluator for a given vertical coordinate.
    # du = f!(u0) at the zero IC IS the surfaced field (constant-RHS D-equation).
    function _lev_eval(lev_coord; F_3d = _make_F3d())
        raw = JSON3.read(read(_LEVMIN_FIXTURE, String))
        ca = Dict("F_3d" => F_3d, "lev_coord" => Float64.(lev_coord))
        ics = Dict("F_surf[$a,$b]" => 0.0 for a in 1:2, b in 1:2)
        f!, u0, p, _, vmap = ESS.build_evaluator(
            raw; model_name = "LevMinSurfaceReduce",
            initial_conditions = ics, const_arrays = ca,
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        return [du[vmap["F_surf[$a,$b]"]] for a in 1:2, b in 1:2]
    end

    # Analytic oracle (no engine): the surface is the F_3d slice at the level
    # whose coordinate is the unique minimum.
    function _surface_oracle(lev_coord; F_3d = _make_F3d())
        lmin = argmin(lev_coord)
        return F_3d[:, :, lmin]
    end
end

@testitem "lev=min 3D->surface reduction end-to-end (esd-khv)" setup = [LevMinReduceSetup] tags = [:reduction, :lev_min, :conformance] begin

    @testset "fixture loads (schema + structural)" begin
        @test isfile(_LEVMIN_FIXTURE)
        @test (ESS.load(_LEVMIN_FIXTURE); true)
    end

    # ACCEPTANCE — output equals the lev=min slice. The canonical fixture places
    # the minimum at an INTERIOR level (lev 2), so the engine must genuinely search
    # for the argmin: a first-slice (lev 1) or last-slice (lev 3) shortcut yields a
    # DIFFERENT field. du = f!(u0) reproduces the analytic slice exactly.
    @testset "surface == lev=min slice (interior minimum, argmin load-bearing)" begin
        lev_coord = [3.0, 1.0, 2.0]                 # unique min at lev 2
        got = _lev_eval(lev_coord)
        @test got == _surface_oracle(lev_coord)     # byte-exact (0-ULP gather)
        @test got == Float64[112 122; 212 222]      # the lev=2 slice (last digit 2)
        # the argmin is load-bearing: neither end-slice equals it
        @test got != _make_F3d()[:, :, 1]           # not the first level
        @test got != _make_F3d()[:, :, 3]           # not the last level
    end

    # The reduction selects value:min WHEREVER the minimum lands — proving it is a
    # genuine argmin reduction, not a hard-coded slice index. First / interior /
    # last, plus realistic monotone pressure coordinates (ascending and descending).
    @testset "argmin generality: minimum at any level / monotone coordinates" begin
        cases = [
            ([1.0, 2.0, 3.0], 1, "ascending ordinal -> first slice (the common surface case)"),
            ([3.0, 1.0, 2.0], 2, "permuted -> interior slice"),
            ([3.0, 2.0, 1.0], 3, "descending -> last slice"),
            ([1000.0, 850.0, 500.0], 3, "descending pressures, min at lev 3"),
            ([500.0, 850.0, 1000.0], 1, "ascending pressures, min at lev 1"),
        ]
        for (lev_coord, lmin, label) in cases
            got = _lev_eval(lev_coord)
            @test got == _make_F3d()[:, :, lmin]            # exact slice at argmin
            @test got == _surface_oracle(lev_coord)         # matches analytic oracle
        end
    end

    # EXACTNESS — `min` returns a lev_coord element bit-for-bit, so the `==`
    # indicator fires at exactly the argmin level and the surfaced value is a
    # verbatim F_3d gather. Holds even for coordinates whose differences are not
    # exactly representable (the compare is element-vs-its-own-min, never an
    # arithmetic result).
    @testset "byte-exact selection for non-round coordinates" begin
        lev_coord = [0.1, 0.3, 0.2]                  # min 0.1 at lev 1; 0.1+0.2 != 0.3 in f64
        F_3d = Float64[a + 0.5b + 0.25l for a in 1:2, b in 1:2, l in 1:3]
        got = _lev_eval(lev_coord; F_3d = F_3d)
        @test got == F_3d[:, :, 1]                   # exact, no tolerance needed
    end
end
