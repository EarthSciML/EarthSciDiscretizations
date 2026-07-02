#!/usr/bin/env julia
# Regenerate the golden for the conservative overlap regridder conformance harness.
#
# The reference-binding golden is the Julia/GeometryOps overlap assembly: build the
# build-once overlap matrix A_ij = polygon_area(intersect_polygon(src_i, tgt_j,
# "spherical"), "spherical"); A_j = column-sums (= dst_areas); F_tgt[j] =
# (1/A_j[j])·Σ_i A_ij·F_src[i]. Per CONFORMANCE_SPEC §5.8 these per-pair areas are
# floating-point clip results (tolerance-based across bindings, NOT byte-identical);
# the PRIMARY gate is the invariants — conservation (Σ_j A_j·F_tgt = Σ_i A_i·F_src)
# and partition-of-unity (Σ_i W_ij = 1) — which the golden also records and which
# the consumer (test/test_regridding_conformance.jl) re-derives THROUGH THE ESS
# ENGINE and asserts.
#
# Needs ESS (intersect_polygon / polygon_area) + GeometryOps (the spherical clip
# backend, a weakdep activated by loading GeometryOps + GeoInterface). Activates a
# throwaway env so the repo deps are untouched.
#
# Run:  julia tests/conformance/regridding/conservative/regenerate_golden.jl

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.add(
        PackageSpec(
            url = "https://github.com/EarthSciML/EarthSciSerialization.git",
            subdir = "packages/EarthSciSerialization.jl", rev = "main",
        ); io = devnull,
    )
    Pkg.add(["GeometryOps", "GeoInterface", "JSON"]; io = devnull)
end
import EarthSciSerialization as ESS
import GeometryOps, GeoInterface   # activates EarthSciSerializationGeometryOpsExt
using JSON

const HERE = @__DIR__
const FIXTURES = joinpath(HERE, "fixtures.json")
const GOLDEN = joinpath(HERE, "golden.json")

ring_matrix(ring) = Float64[ring[k][d] for k in 1:length(ring), d in 1:2]

function main()
    spec = JSON.parsefile(FIXTURES)
    src = [ring_matrix(r) for r in spec["src_polys"]]
    tgt = [ring_matrix(r) for r in spec["tgt_polys"]]
    F_src = Float64.(spec["F_src"])
    nsrc = length(src); ntgt = length(tgt)

    A_ij = zeros(Float64, nsrc, ntgt)
    for i in 1:nsrc, j in 1:ntgt
        ring = ESS.intersect_polygon(src[i], tgt[j], "spherical")
        size(ring, 1) < 3 && continue
        A_ij[i, j] = ESS.polygon_area(ring, "spherical")
    end
    A_j = vec(sum(A_ij; dims = 1))        # column-sums = dst_areas
    A_i = vec(sum(A_ij; dims = 2))        # row-sums = per-source overlap mass
    F_tgt = [sum(A_ij[i, j] * F_src[i] for i in 1:nsrc) / A_j[j] for j in 1:ntgt]

    src_mass = sum(A_i .* F_src)
    tgt_mass = sum(A_j .* F_tgt)
    pou = [sum(A_ij[i, j] for i in 1:nsrc) / A_j[j] for j in 1:ntgt]

    out = Dict{String, Any}(
        "model" => "conservative",
        "reference_binding" => "julia",
        "manifold" => "spherical",
        "A_j" => A_j,
        "F_tgt" => F_tgt,
        "partition_of_unity" => pou,        # == 1 per target, exact by construction
        "source_mass" => src_mass,          # Σ_i A_i·F_src
        "target_mass" => tgt_mass,          # Σ_j A_j·F_tgt (== source_mass)
    )
    open(GOLDEN, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote $GOLDEN")
    println("  A_j=$A_j")
    println("  F_tgt=$F_tgt")
    println("  source_mass=$src_mass  target_mass=$tgt_mass  (Δ=$(abs(src_mass - tgt_mass)))")
    println("  partition_of_unity=$pou")
    return
end

main()
