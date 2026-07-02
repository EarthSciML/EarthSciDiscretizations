# Standalone verification (plain Test) that conservative_regrid_overlap_join_computed.esm
# COMPUTES the overlap matrix A_ij internally — given only the source/target polygon
# geometry, with NO host-supplied A_ij/dst_areas — and reproduces the host-A_ij
# conformance fixture's A_j/F_tgt on the IDENTICAL 3×3 spherical case.
#
# Run through the ESS pde_sim_adapter env (which pins ESS + GeometryOps + GeoInterface
# + JSON3):  cd <ESS>/packages/EarthSciSerialization.jl/scripts/pde_sim_adapter &&
#            julia --project=. -e 'include(".../conservative_regrid_overlap_join_computed_verify.jl")'
using Test
import EarthSciSerialization as ESS
import GeometryOps, GeoInterface          # activate the spherical clip backend extension
const JSON3 = ESS.JSON3

const FIX = joinpath(@__DIR__, "conservative_regrid_overlap_join.esm")
const CMP = joinpath(@__DIR__, "conservative_regrid_overlap_join_computed.esm")

# The fixture's exact 3×3 spherical case (conservative_regrid_overlap_join_conformance_test.jl).
_rect(x0, x1, y0, y1) = [x0 y0; x1 y0; x1 y1; x0 y1]      # CCW lon/lat ring (degrees)
const _SRC = [_rect(0, 1, 0, 1), _rect(1, 2, 0, 1), _rect(2, 3, 0, 1)]
const _TGT = [_rect(0, 1.5, 0, 1), _rect(1.5, 2, 0, 1), _rect(2, 3, 0, 1)]
const _SRC_LON = [0.0, 1.0, 2.0]; const _SRC_LAT = [0.0, 0.0, 0.0]
const _TGT_LON = [0.0, 1.5, 2.0]; const _TGT_LAT = [0.0, 0.0, 0.0]
const _F_SRC = [10.0, 20.0, 30.0]
const _DX = 2.0; const _DY = 2.0

# Oracle: build A_ij the way the fixture's host does — the spherical clip + VOS area.
function oracle()
    A = zeros(3, 3)
    for i in 1:3, j in 1:3
        ring = ESS.intersect_polygon(_SRC[i], _TGT[j], "spherical")
        size(ring, 1) < 3 && continue
        A[i, j] = ESS.polygon_area(ring, "spherical")
    end
    dst = vec(sum(A; dims = 1))
    F = [sum(A[i, j] * _F_SRC[i] for i in 1:3) / dst[j] for j in 1:3]
    return A, dst, F
end

# Pack the 3 cell rings into a [cells,4,2] polygon array for the computed component.
polys(rings) = (P = zeros(length(rings), 4, 2);
                for (k, r) in enumerate(rings); P[k, :, :] = r; end; P)

@testset "computed-A_ij conservative regrid reproduces the host-A_ij fixture" begin
    A_ij_oracle, dst_oracle, F_oracle = oracle()

    raw = JSON3.read(read(CMP, String))
    ca = Dict(
        "src_poly" => polys(_SRC), "tgt_poly" => polys(_TGT), "F_src" => _F_SRC,
        "src_lon" => _SRC_LON, "src_lat" => _SRC_LAT,
        "tgt_lon" => _TGT_LON, "tgt_lat" => _TGT_LAT,
    )                                                  # NOTE: no A_ij, no dst_areas
    ics = Dict("A_j[1]" => 0.0, "A_j[2]" => 0.0, "A_j[3]" => 0.0,
               "F_tgt[1]" => 0.0, "F_tgt[2]" => 0.0, "F_tgt[3]" => 0.0)
    f!, u0, p, _t, vmap = ESS.build_evaluator(
        raw; model_name = "ConservativeRegridOverlapJoinComputed",
        initial_conditions = ics, const_arrays = ca,
        parameter_overrides = Dict("dx" => _DX, "dy" => _DY, "atol" => 1.0e-15))
    du = similar(u0); f!(du, u0, p, 0.0)
    A_j = [du[vmap["A_j[$j]"]] for j in 1:3]
    F_tgt = [du[vmap["F_tgt[$j]"]] for j in 1:3]

    println("  host-fixture oracle  A_j  = ", round.(dst_oracle; digits = 8))
    println("  COMPUTED component   A_j  = ", round.(A_j; digits = 8))
    println("  host-fixture oracle  F_tgt= ", round.(F_oracle; digits = 6))
    println("  COMPUTED component   F_tgt= ", round.(F_tgt; digits = 6))

    @test isapprox(A_j, dst_oracle; rtol = 1e-9, atol = 1e-14)      # computed denominator = oracle
    @test isapprox(F_tgt, F_oracle; rtol = 1e-9, atol = 1e-12)      # full conservative regrid = oracle
    # conservation: Σ_j A_j·F_tgt = Σ_i A_i·F_src (overlap-mass row-sums)
    src_areas = vec(sum(A_ij_oracle; dims = 2))
    @test isapprox(sum(A_j .* F_tgt), sum(src_areas .* _F_SRC); rtol = 1e-9)
    # partition of unity: Σ_i A_ij / A_j = 1 per target cell
    for j in 1:3
        @test isapprox(sum(A_ij_oracle[i, j] for i in 1:3) / A_j[j], 1.0; rtol = 1e-9)
    end
    println("  ✓ A_ij computed from polygon geometry (no host overlap matrix); ",
            "reproduces the fixture's conservative F_tgt, conservative + partition-of-unity.")
end
