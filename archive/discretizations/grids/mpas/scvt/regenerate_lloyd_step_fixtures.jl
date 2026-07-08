#!/usr/bin/env julia
# Regenerate the canonical MPAS-SCVT one-iteration Lloyd STEP `.esm` golden
# (bead esd-e5m.2 / D2) from the Julia reference binding.
#
# Per AGENTS.md ("golden regeneration drives the canonical pipeline"), every
# value is produced THROUGH the single pathway: the background quadrature point
# set + density-weighted measure come from the D1 FAQ (eval_coeff of the same
# centroid / rho / mass ASTs as background_quadrature.esm), and the SCVT STEP
# itself — assign (E1 argmin), the grouped density-weighted centroid num/den (E2
# grouped sum_product), and the elementwise move centroid = num/den — is the
# end-to-end output of `EarthSciSerialization.materialize_value_invention`
# driving lloyd_step.esm. No shadow evaluator, no binding-local shortcut.
#
# Canonical instance: the level-0 composition of D1 (n_cells = 20 background
# quadrature points = the level-0 DUO primal cells) with the icosahedral
# generator set (n_generators = 12 = the level-0 icosahedron vertices). The 20
# face centroids are nearest, by the §5.7 smallest-id tie-break, to the smaller
# of their 3 corner generators, so 4 of the 12 generators are unattended at this
# coarse quadrature (den = 0, centroid = NaN — the empty-group identity); a
# finer background attends every generator (see the convergence test). assign is
# density-independent (a pure distance argmin), so it is shared by both density
# samples; den / num / centroid are pinned for each of:
#   rho == 1            — the uniform-density (CVT) regression.
#   rho == 2 + z        — a smooth latitude-graded density (z = unit-centroid
#                         z-component; rho in [1,3]) driving variable resolution.
#
# Run from the repo root:
#     julia --project=. discretizations/grids/mpas/scvt/regenerate_lloyd_step_fixtures.jl

using EarthSciDiscretizations
const ESD = EarthSciDiscretizations
using EarthSciSerialization
const ESS = EarthSciSerialization
using JSON

const REPO_ROOT = dirname(dirname(dirname(dirname(dirname(@__FILE__)))))
const SCVT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt")
const OUT_DIR = joinpath(SCVT_DIR, "fixtures", "canonical")

mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

# Centroid AST over the three corner unit vectors a,b,c (= primal_geometry ux/uy/uz),
# identical to background_quadrature.esm's bg_coord (D1) — so D2's quadrature
# point set is byte-identical to D1's.
const MX = mk("+", "a1", "b1", "c1")
const MY = mk("+", "a2", "b2", "c2")
const MZ = mk("+", "a3", "b3", "c3")
const NRM = mk("sqrt", mk("+", mk("*", MX, MX), mk("*", MY, MY), mk("*", MZ, MZ)))
centroid_comp(m) = mk("/", m, NRM)
const UCOMP = (centroid_comp(MX), centroid_comp(MY), centroid_comp(MZ))

const RHO_AST = mk("+", 2.0, "z")          # rho(x,y,z) = 2 + z
const MASS_AST = mk("*", "rho", "area")    # bg_mass = rho * area (D1 q-1)

corner_binding(V, i, j, k) = Dict{String, Float64}(
    "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
    "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
    "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k],
)

# NaN centroid (an unattended generator) serialises as JSON null — the golden
# stays valid JSON; the test maps null back to an isnan assertion.
nan_to_null(v) = Any[isnan(x) ? nothing : x for x in v]

function canonical_bytes(d::AbstractDict)
    io = IOBuffer()
    JSON.print(io, d, 2)
    print(io, "\n")
    return String(take!(io))
end

function build_level0()
    level = 0
    g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
    V, F = ESD.duo_subdivide_faq(Float64, level)
    R = g.R
    Nc = ESD.n_cells(g)               # 20 background quadrature points
    Ng = size(V, 2)                   # 12 icosahedral-vertex generators

    # Background quadrature points (unit-sphere face centroids) + spherical-tri areas.
    bg_coord = Matrix{Float64}(undef, Nc, 3)
    area = Vector{Float64}(undef, Nc)
    for c in 1:Nc
        b = corner_binding(V, F[1, c], F[2, c], F[3, c])
        for d in 1:3
            bg_coord[c, d] = ESD.eval_coeff(UCOMP[d], b)
        end
        area[c] = g.area[c]
    end
    # Generators = the 12 unit icosahedral vertices.
    gen = Matrix{Float64}(undef, Ng, 3)
    for v in 1:Ng, d in 1:3
        gen[v, d] = V[d, v]
    end

    # Two density-weighted measures bg_mass = rho * area, sampled through eval_coeff.
    rho_uniform = ones(Float64, Nc)
    rho_sampled = [ESD.eval_coeff(RHO_AST, Dict{String, Float64}("z" => bg_coord[c, 3])) for c in 1:Nc]
    bg_mass(rho) = [ESD.eval_coeff(MASS_AST, Dict{String, Float64}("rho" => rho[c], "area" => area[c])) for c in 1:Nc]

    # Drive the SCVT step end-to-end through the value-invention front-door.
    raw = JSON.parsefile(joinpath(SCVT_DIR, "lloyd_step.esm"); dicttype = Dict{String, Any})
    mj = ESS._select_model_json(raw, "ScvtLloydStep")
    @assert mj["index_sets"]["cells"]["size"] == Nc
    @assert mj["index_sets"]["generators"]["size"] == Ng

    function step(mass)
        ca = Dict("bg_coord" => bg_coord, "bg_mass" => mass, "gen" => gen)
        vi = ESS.materialize_value_invention(mj, ca, Dict{String, Float64}())
        return vi
    end

    mass_u = bg_mass(rho_uniform)
    mass_s = bg_mass(rho_sampled)
    vi_u = step(mass_u)
    vi_s = step(mass_s)
    # assign is a pure distance argmin — density-independent (shared).
    @assert vi_u.assignments["assign"] == vi_s.assignments["assign"]

    density_block(mass, vi) = Dict{String, Any}(
        "bg_mass" => mass,
        "den" => vi.groups["den"],
        "num_x" => vi.groups["num_x"],
        "num_y" => vi.groups["num_y"],
        "num_z" => vi.groups["num_z"],
        "centroid_x" => nan_to_null(vi.groups["centroid_x"]),
        "centroid_y" => nan_to_null(vi.groups["centroid_y"]),
        "centroid_z" => nan_to_null(vi.groups["centroid_z"]),
    )

    return Dict{String, Any}(
        "description" => "Canonical cross-binding contract (bead esd-e5m.2 / D2) for the MPAS-SCVT one-iteration Lloyd STEP, the level-0 composition of D1 (n_cells = 20 background quadrature points) with the icosahedral generator set (n_generators = 12). Produced end-to-end through EarthSciSerialization.materialize_value_invention driving lloyd_step.esm. CONTRACT: the integer `assign` buffer (E1 argmin, §5.7 rule-6 smallest-id tie-break) is BYTE-identical across bindings; the grouped/derived float buffers (den, num_{x,y,z}, centroid_{x,y,z}) ride the §5.8 tolerance contract across bindings and are reproduced bit-for-bit by the Julia reference. assign is density-independent (a pure distance argmin) so it is shared; den/num/centroid are pinned for the uniform density (rho==1) and a sampled latitude-graded density (rho = 2 + bg_coord_z). A centroid of null marks an unattended generator (no assigned point at this coarse quadrature: den=0, centroid=NaN=0/0). Pinned by test/test_mpas_scvt_lloyd_step_faq.jl. Regenerate with discretizations/grids/mpas/scvt/regenerate_lloyd_step_fixtures.jl.",
        "level" => level,
        "n_cells" => Nc,
        "n_generators" => Ng,
        "R" => R,
        "density_sample" => "rho(x,y,z) = 2 + z  (z = unit-centroid z-component; rho in [1,3], always positive)",
        "bg_coord_x" => bg_coord[:, 1],
        "bg_coord_y" => bg_coord[:, 2],
        "bg_coord_z" => bg_coord[:, 3],
        "gen_x" => gen[:, 1],
        "gen_y" => gen[:, 2],
        "gen_z" => gen[:, 3],
        "area" => area,
        "assign" => vi_u.assignments["assign"],
        "rho_sampled" => rho_sampled,
        "uniform" => density_block(mass_u, vi_u),
        "sampled" => density_block(mass_s, vi_s),
    )
end

function main()
    isdir(OUT_DIR) || mkpath(OUT_DIR)
    d = build_level0()
    path = joinpath(OUT_DIR, "lloyd_step_level0.json")
    write(path, canonical_bytes(d))
    println("wrote $(path) (n_cells=$(d["n_cells"]), n_generators=$(d["n_generators"]))")
    return
end

main()
