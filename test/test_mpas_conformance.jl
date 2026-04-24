@testsnippet MpasConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON

    # This conformance harness lives at tests/conformance/grids/mpas/ and
    # is binding-neutral: fixtures.json carries the full MPAS mesh arrays
    # inline (since NetCDF mesh files are not bundled in this repo). Each
    # binding constructs an MpasMeshData from those arrays, builds a
    # grid, queries accessors at pinned points, and compares to the
    # committed golden.
    #
    # The mpas reference binding is interim Python until bead dsc-7j0
    # (Julia mpas runtime) lands. Once the Julia `mpas` accessor runtime
    # is available at `EarthSciDiscretizations.grids.mpas` (or the flat
    # `EarthSciDiscretizations.mpas`), this file flips from the runtime
    # `@test_skip` guard below to a first-class accessor comparison that
    # mirrors `python/tests/test_mpas_conformance.py`.

    const MPAS_HARNESS_DIR =
        joinpath(@__DIR__, "..", "tests", "conformance", "grids", "mpas")
    const MPAS_JULIA_AVAILABLE = isdefined(EarthSciDiscretizations, :MpasGrid)
end

@testitem "mpas cross-language conformance" setup = [MpasConformanceSetup] tags =
    [:conformance, :grid, :mpas] begin

    # Always assert the harness files are present so structural drift is
    # caught even before the Julia binding lands.
    fixtures_path = joinpath(MPAS_HARNESS_DIR, "fixtures.json")
    @test isfile(fixtures_path)
    spec = JSON.parsefile(fixtures_path)
    @test spec["family"] == "mpas"
    @test haskey(spec, "tolerance")
    @test haskey(spec, "fixtures")

    for fixture in spec["fixtures"]
        name = fixture["name"]
        golden_path = joinpath(MPAS_HARNESS_DIR, "golden", "$name.json")
        @test isfile(golden_path)
    end

    if !MPAS_JULIA_AVAILABLE
        # Julia mpas runtime is still in flight on bead dsc-7j0. Mark the
        # accessor-comparison step as skipped-pending-binding. See
        # tests/conformance/grids/mpas/README.md for the rollout plan.
        @test_skip "julia mpas runtime not yet available (bead dsc-7j0)"
    else
        # When dsc-7j0 lands, flip this branch to a full accessor
        # comparison mirroring python/tests/test_mpas_conformance.py:
        # for each fixture, build an MpasGrid from fixtures["mesh"],
        # probe cell_centers / neighbors / cell_area / edge_length /
        # metric_eval at the pinned query indices, and compare to
        # golden at REL_TOL = spec["tolerance"]["relative"].
        error(
            "mpas conformance: Julia binding present but accessor test body not " *
            "yet implemented. Port python/tests/test_mpas_conformance.py here.",
        )
    end
end
