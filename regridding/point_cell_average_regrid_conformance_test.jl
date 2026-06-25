# Point cell-averaging regridder — Julia evaluator conformance.
#
# ESD owns this declarative regridder PROGRAM (point_cell_average_regrid.esm,
# co-located here); the EarthSciSerialization engine (an ESD dependency)
# EVALUATES it. Sibling of conservative_regrid_overlap_join.esm — same bin-Skolem
# broad-phase idiom (floor -> skolem bin keys -> distinct equi-join), but the
# source is SCATTERED POINT STATIONS, not gridded cells, and the assembly is a
# per-cell AVERAGE with a configurable missing_value fill rather than an
# area-weighted normalize.
#
# RFC pure-io-data-loaders §5.2 (point loaders / OpenAQ); bead esd-47z.4. The
# OpenAQ cell-averaging method (EarthSciData.jl's `_build_cell_station_map` +
# per-cell mean) as ONE document, driven through build_evaluator with NO
# pre-baked bins:
#
#   floor -> skolem station/cell bin keys -> distinct assignment set   (broad phase)
#   -> count group-by (Σ 1)                                            (per-cell tally)
#   -> average group-by (Σ val/count_denom), ifelse missing fill       (assembly)
#
# Each assembled state is a constant-RHS D-equation from a zero IC, so the
# derivative du = f!(u0) IS the assembled value — exact, no integrator (the
# conservative_regrid assembly precedent). The fixture is supplied
# station_val / count_denom / coords as const_arrays, so per CONFORMANCE_SPEC
# §5.8 the other bindings schema-/structurally-validate it while Julia evaluates
# it here.

@testsnippet PointCellAverageSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach it
    # through the dependency so the raw document is the exact type the
    # value-invention front-door expects.
    const JSON3 = ESS.JSON3

    const _PT_FIXTURE = joinpath(@__DIR__, "point_cell_average_regrid.esm")

    # 3 scattered stations over a target grid of 3 unit cells (dx=dy=1) in a thin
    # equatorial band. The cell representative is its lower-left corner, so
    # floor(cell_lon)=column index. With dx=1 the bins are:
    #   station 1 (0.3,0.5) -> bin (0,0)   station 2 (0.7,0.2) -> bin (0,0)
    #   station 3 (1.5,0.5) -> bin (1,0)
    #   cell 1 (0,0) -> bin (0,0)  cell 2 (1,0) -> bin (1,0)  cell 3 (2,0) -> bin (2,0)
    # so cell 1 gets stations {1,2} (mean (10+20)/2 = 15), cell 2 gets {3}
    # (mean 30), cell 3 gets NONE -> missing_value. The broad phase is
    # LOAD-BEARING: the assignment is the within-bin block, not the full 3x3.
    const _ST_LON = [0.3, 0.7, 1.5]; const _ST_LAT = [0.5, 0.2, 0.5]
    const _CELL_LON = [0.0, 1.0, 2.0]; const _CELL_LAT = [0.0, 0.0, 0.0]
    const _ST_VAL = [10.0, 20.0, 30.0]
    const _DX = 1.0; const _DY = 1.0
    const _MISSING = -999.0
    # Build-once per-cell contributing-station count (EarthSciData's station tally)
    # and the expected per-cell mean — the host oracle the FAQs must reproduce.
    const _COUNT_DENOM = [2.0, 1.0, 0.0]
    const _F_EXPECT = [15.0, 30.0, _MISSING]

    _pt_const_arrays() = Dict(
        "station_lon" => _ST_LON, "station_lat" => _ST_LAT,
        "cell_lon" => _CELL_LON, "cell_lat" => _CELL_LAT,
        "station_val" => _ST_VAL, "count_denom" => _COUNT_DENOM,
    )

    # Drive the point-regridder fixture through build_evaluator, returning
    # (count, F_tgt). du = f!(u0) at the zero IC IS the assembled value for each
    # constant-RHS D-equation.
    function _pt_eval(; count_denom = _COUNT_DENOM, missing_value = _MISSING)
        # The VALUE-INVENTION broad phase (station_bin/cell_bin/assignments) is
        # materialised only on the AbstractDict front-door, so drive the RAW
        # document through build_evaluator (NOT the typed EsmFile path, which
        # skips it).
        raw = JSON3.read(read(_PT_FIXTURE, String))
        ca = merge(_pt_const_arrays(), Dict("count_denom" => count_denom))
        ics = Dict(
            ("count[$j]" => 0.0 for j in 1:3)...,
            ("F_tgt[$j]" => 0.0 for j in 1:3)...,
        )
        f!, u0, p, _, vmap = ESS.build_evaluator(
            raw; model_name = "PointCellAverageRegrid",
            initial_conditions = ics, const_arrays = ca,
            parameter_overrides = Dict("dx" => _DX, "dy" => _DY,
                                       "missing_value" => missing_value),
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        count = [du[vmap["count[$j]"]] for j in 1:3]
        F_tgt = [du[vmap["F_tgt[$j]"]] for j in 1:3]
        return count, F_tgt
    end
end

@testitem "point cell-averaging regridder end-to-end (esd-47z.4)" setup = [PointCellAverageSetup] tags = [:regridding, :point, :conformance] begin

    @testset "fixture loads (schema + structural)" begin
        @test isfile(_PT_FIXTURE)
        @test (ESS.load(_PT_FIXTURE); true)
    end

    # (1) BROAD PHASE — the bin-Skolem assignment set is materialised through the
    # value-invention front-door from the COMPUTED bins (no pre-baked members).
    # dx=1 => stations {1,2} share cell-1's bin, station 3 is in cell-2's bin,
    # cell 3 gets nothing.
    @testset "broad phase: assignment set = computed bin-Skolem equi-join" begin
        mj = ESS._select_model_json(
            JSON3.read(read(_PT_FIXTURE, String)),
            "PointCellAverageRegrid",
        )
        vi = ESS.materialize_value_invention(
            mj, Dict(
                "station_lon" => _ST_LON, "station_lat" => _ST_LAT,
                "cell_lon" => _CELL_LON, "cell_lat" => _CELL_LAT,
            ),
            Dict("dx" => _DX, "dy" => _DY, "missing_value" => _MISSING),
        )
        # surviving (station i, cell j) assignments in §5.7 sorted order
        @test vi.members["assignment_set"] == [(1, 1), (2, 1), (3, 2)]
        @test vi.extents["assignment_set"] == 3
        @test vi.vi_var_names == Set(["station_bin", "cell_bin", "pair_exists"])
        # the materialised bin buffers the F3 join gates on
        @test vi.maps["station_bin"][1] == vi.maps["station_bin"][2]   # stations 1,2 same bin
        @test vi.maps["station_bin"][1] != vi.maps["station_bin"][3]   # station 3 different
        @test vi.maps["station_bin"][3] == vi.maps["cell_bin"][2]      # station 3 in cell 2
        @test vi.map_sets["station_bin"] == "stations"
        @test vi.map_sets["cell_bin"] == "cells"
    end

    # (2) COUNT group-by over the BUFFER-GATED assignment join reproduces the
    # build-once count_denom (the F3 join.on [[station_bin,cell_bin]] capability).
    @testset "count reproduces count_denom (buffer-gated join)" begin
        count, _ = _pt_eval()
        @test count ≈ _COUNT_DENOM
        @test count == [2.0, 1.0, 0.0]
    end

    # ACCEPTANCE 1 — CELL MEAN of contributing stations. cell 1 = mean(10,20) =
    # 15, cell 2 = mean(30) = 30.
    @testset "average: cell mean of contributing stations" begin
        _, F_tgt = _pt_eval()
        @test F_tgt[1] ≈ 15.0
        @test F_tgt[2] ≈ 30.0
    end

    # ACCEPTANCE 2 — configurable missing_value for a cell with NO station. cell 3
    # is empty -> the configured fill, NOT 0 and NOT NaN.
    @testset "missing_value: empty cell returns the configured fill" begin
        _, F_tgt = _pt_eval()
        @test F_tgt[3] == _MISSING
        @test F_tgt ≈ _F_EXPECT
        # the fill is genuinely configurable, not a hard-coded sentinel
        _, F_alt = _pt_eval(missing_value = 1.0e20)
        @test F_alt[3] == 1.0e20
        @test F_alt[1] ≈ 15.0   # non-empty cells unaffected by the fill value
    end

    # ACCEPTANCE 3 — PARTITION-OF-UNITY: the average is a weighted sum with
    # uniform weights W_ij = 1/count_denom[j] over the contributing stations, so
    # Σ_i W_ij = 1 exactly (the mean is a partition of unity). Equivalent to:
    # the regridded value equals the arithmetic mean of the assigned values.
    @testset "PARTITION-OF-UNITY: Σ_i W_ij = 1 over contributing stations" begin
        _, F_tgt = _pt_eval()
        # cell 1: (1/2 + 1/2) = 1 -> F = (10+20)/2; cell 2: (1/1) = 1 -> F = 30
        @test F_tgt[1] ≈ sum(_ST_VAL[1:2]) / _COUNT_DENOM[1]
        @test F_tgt[2] ≈ _ST_VAL[3] / _COUNT_DENOM[2]
    end

    # The BUFFER-GATED join is load-bearing: a station whose bin does NOT match a
    # cell must NOT contribute. Corrupting count_denom for the empty cell to a
    # positive value must STILL yield no contribution (Σ over an empty assignment
    # block = 0), proving the broad phase — not the denominator — selects terms.
    @testset "buffer-gated join: empty cell has no contributing terms" begin
        # count_denom[3] = 1 (a lie: cell 3 has no station). The average body
        # Σ_i station_val[i]/count_denom[3] still ranges an EMPTY assignment block
        # -> 0, and count[3] (the honest tally) stays 0.
        count, F_tgt = _pt_eval(count_denom = [2.0, 1.0, 1.0])
        @test count[3] == 0.0                 # honest group-by tally: still empty
        @test F_tgt[3] ≈ 0.0                  # empty Σ, divided by the lie -> 0
        @test F_tgt[1] ≈ 15.0                 # real cells unchanged
        @test F_tgt[2] ≈ 30.0
    end
end
