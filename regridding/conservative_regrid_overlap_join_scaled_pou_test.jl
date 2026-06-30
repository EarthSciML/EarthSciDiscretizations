# Scaled partition-of-unity regression for conservative_regrid_overlap_join_computed.esm.
#
# The 3×3 conformance fixture (conservative_regrid_overlap_join_computed_verify.jl) has
# no EDGE-ADJACENT cell pairs, so it never exercised the failure this guards: at scale a
# spherical clip emits a ~1e-9 sliver for cells that share only an EDGE (zero true
# overlap). The build-once denominator A_j_w must contract over the SAME broad-phase
# bin-skolem candidate set + sub-atol sliver filter as the apply numerator; a DENSE
# Σ_i A_ij keeps those slivers and biases F_tgt by O(sliver/A_j), breaking the partition
# of unity. (Joint fix: the component declares A_j_w's join+filter, and ESS threads the
# value-invention bin buffers + scalar param overrides into the setup-time geometry env so
# the gate resolves — without it the join silently no-ops and the filter throws on `atol`.)
#
# A 2-D NESTED grid: each target cell falls fully inside exactly one source cell, and the
# refinement places a band of target cells with an edge exactly on every interior source
# boundary — the cells the dense denominator used to bias. The conservative answer is
# F_tgt = the containing source's F_src, exactly.
#
# Run through the ESS pde_sim_adapter env:
#   cd <ESS>/packages/EarthSciSerialization.jl/scripts/pde_sim_adapter &&
#   julia --project=. .../conservative_regrid_overlap_join_scaled_pou_test.jl
using Test
import EarthSciSerialization as ESS
import GeometryOps, GeoInterface
const JSON3 = ESS.JSON3
const CMP = joinpath(@__DIR__, "conservative_regrid_overlap_join_computed.esm")

_rect(x0, x1, y0, y1) = [x0 y0; x1 y0; x1 y1; x0 y1]
polys(rings) = (P = zeros(length(rings), 4, 2);
                for (k, r) in enumerate(rings); P[k, :, :] = r; end; P)

# source SXY×SXY cells of size HS; target (SXY·TR)×(SXY·TR) cells of size HS/TR.
function build(SXY, TR)
    HS = 1.0; HT = HS / TR; TXY = SXY * TR; NS = SXY^2; NT = TXY^2
    src = [_rect((c-1)*HS, c*HS, (r-1)*HS, r*HS) for r in 1:SXY for c in 1:SXY]
    tgt = [_rect((c-1)*HT, c*HT, (r-1)*HT, r*HT) for r in 1:TXY for c in 1:TXY]
    src_lon = Float64[((k-1)%SXY)*HS for k in 1:NS]; src_lat = Float64[((k-1)÷SXY)*HS for k in 1:NS]
    tgt_lon = Float64[((k-1)%TXY)*HT for k in 1:NT]; tgt_lat = Float64[((k-1)÷TXY)*HT for k in 1:NT]
    F_src = Float64[10k for k in 1:NS]
    (; HS, HT, TXY, NS, NT, src, tgt, src_lon, src_lat, tgt_lon, tgt_lat, F_src)
end

# Analytically-correct reference: each target nests fully in exactly ONE source, so the
# conservative value is that source's F_src exactly and the area is the target's own
# (spherical) cell area. NOTE: a full-dense Σ-over-ALL-pairs oracle would itself be
# BIASED by the ~1e-9 edge-adjacent slivers (it is the very computation the broad-phase
# denominator removes), so the correct reference is the container value, not the dense sum.
function oracle(g)
    SXY = isqrt(g.NS)
    dst = zeros(g.NT); F = zeros(g.NT)
    for j in 1:g.NT
        cx = Int(floor(g.tgt_lon[j] / g.HS)); cy = Int(floor(g.tgt_lat[j] / g.HS))
        i = cy * SXY + cx + 1                                    # containing source (col,row)
        ring = ESS.intersect_polygon(g.src[i], g.tgt[j], "spherical")
        size(ring, 1) >= 3 && (dst[j] = ESS.polygon_area(ring, "spherical"))
        F[j] = g.F_src[i]
    end
    dst, F
end

function run_case(SXY, TR)
    g = build(SXY, TR)
    dst_or, F_or = oracle(g)
    raw = JSON3.read(read(CMP, String), Dict{String,Any})
    model = raw["models"]["ConservativeRegridOverlapJoinComputed"]
    model["index_sets"]["src_cells"]["size"] = g.NS
    model["index_sets"]["tgt_cells"]["size"] = g.NT
    for eq in model["equations"]                       # widen the [1,3] ODE LHS ranges
        lhs = eq["lhs"]
        lhs isa AbstractDict && get(lhs, "op", "") == "aggregate" &&
            haskey(lhs["ranges"], "j") && (lhs["ranges"]["j"] = Any[1, g.NT])
    end
    ca = Dict("src_poly"=>polys(g.src), "tgt_poly"=>polys(g.tgt), "F_src"=>g.F_src,
              "src_lon"=>g.src_lon, "src_lat"=>g.src_lat, "tgt_lon"=>g.tgt_lon, "tgt_lat"=>g.tgt_lat)
    ics = Dict{String,Float64}()
    for j in 1:g.NT; ics["A_j[$j]"]=0.0; ics["F_tgt[$j]"]=0.0; end
    f!, u0, p, _t, vmap = ESS.build_evaluator(raw;
        model_name="ConservativeRegridOverlapJoinComputed",
        initial_conditions=ics, const_arrays=ca,
        parameter_overrides=Dict("dx"=>g.HS, "dy"=>g.HS, "atol"=>1.0e-15))
    du = similar(u0); f!(du, u0, p, 0.0)
    A_j   = [du[vmap["A_j[$j]"]] for j in 1:g.NT]
    F_tgt = [du[vmap["F_tgt[$j]"]] for j in 1:g.NT]
    (; g, dst_or, F_or, A_j, F_tgt)
end

@testset "scaled conservative regrid: edge-adjacent partition of unity" begin
    for (SXY, TR) in ((2, 3), (4, 5))
        r = run_case(SXY, TR)
        nedge = count(j -> any(isapprox(r.g.tgt_lon[j], b; atol=1e-9) ||
                               isapprox(r.g.tgt_lat[j], b; atol=1e-9)
                               for b in (1.0:1.0:(SXY-1))), 1:r.g.NT)
        @test isapprox(r.A_j, r.dst_or; rtol=1e-9, atol=1e-12)        # broad-phase area = oracle
        @test isapprox(r.F_tgt, r.F_or; rtol=1e-7, atol=1e-9)         # exact conservative value
        # partition of unity per target: Σ_i A_ij / A_j_w == 1 (here A_j == A_j_w).
        @test all(j -> r.A_j[j] == 0 || isapprox(r.F_tgt[j], r.F_or[j]; rtol=1e-7), 1:r.g.NT)
        println("  $(SXY)×$(SXY) → $(r.g.TXY)×$(r.g.TXY) ($(r.g.NS)→$(r.g.NT) cells, ",
                "$nedge edge-adjacent): max |F_tgt−oracle| = ",
                round(maximum(abs, r.F_tgt .- r.F_or); sigdigits=3))
    end
    println("  ✓ A_j_w broad-phase denominator keeps the partition of unity exact at scale ",
            "(a dense Σ_i biased ~", "1e-4 on every edge-adjacent target).")
end
