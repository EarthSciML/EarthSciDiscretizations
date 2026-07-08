#!/usr/bin/env julia
# Regenerate the canonical MPAS-SCVT background-quadrature `.esm` FAQ goldens
# (bead esd-e5m.1 / D1) from the Julia reference binding.
#
# Per AGENTS.md ("golden regeneration drives the canonical pipeline"), every
# value is produced THROUGH the single-pathway evaluator
# `EarthSciDiscretizations.eval_coeff` (a passthrough to the ESS engine) — never
# a binding-local shortcut. Specifically:
#   - bg_coord (the quadrature point) = the unit-sphere primal-cell centroid,
#     eval_coeff of the centroid AST (= primal_geometry.esm cell_cart / R),
#   - rho             = eval_coeff of the density AST sampled at bg_coord,
#   - bg_mass         = eval_coeff(rho * area)            (the FAQ q-1 product),
#   - bg_moment       = eval_coeff(bg_mass * bg_coord)    (the FAQ q-2 product).
# `area` is the canonical DUO primal-cell area from build_duo_grid (itself proven
# == primal_geometry.esm to 0 ULP by test/test_duo_primal_geometry_faq.jl).
#
# Two density samples are pinned per level:
#   rho == 1            — the uniform-density regression (SCVT -> CVT).
#   rho == 2 + z        — a smooth latitude-graded density (z = unit-centroid
#                         z-component; rho in [1,3], always positive), the
#                         variable-resolution sampling demonstration.
#
# Run from the repo root:
#     julia --project=. discretizations/grids/mpas/scvt/regenerate_fixtures.jl

using EarthSciDiscretizations
const ESD = EarthSciDiscretizations
using JSON

const REPO_ROOT = dirname(dirname(dirname(dirname(dirname(@__FILE__)))))
const OUT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt", "fixtures", "canonical")

mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))

# Centroid AST over the three corner unit vectors a,b,c (= primal_geometry ux/uy/uz).
const MX = mk("+", "a1", "b1", "c1")
const MY = mk("+", "a2", "b2", "c2")
const MZ = mk("+", "a3", "b3", "c3")
const NRM = mk("sqrt", mk("+", mk("*", MX, MX), mk("*", MY, MY), mk("*", MZ, MZ)))
centroid_comp(m) = mk("/", m, NRM)
const UCOMP = (centroid_comp(MX), centroid_comp(MY), centroid_comp(MZ))

# Density sampling AST: rho(x,y,z) = 2 + z (a function of the cartesian point).
const RHO_AST = mk("+", 2.0, "z")

# FAQ integrand product ASTs (background_quadrature.esm q-1 / q-2).
const MASS_AST = mk("*", "rho", "area")          # bg_mass   = rho * area
const MOMENT_AST = mk("*", "mass", "coord")      # bg_moment = bg_mass * bg_coord

corner_binding(V, i, j, k) = Dict{String, Float64}(
    "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
    "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
    "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k],
)

function canonical_bytes(d::AbstractDict)
    io = IOBuffer()
    JSON.print(io, d, 2)
    print(io, "\n")
    return String(take!(io))
end

function build_level(level::Int)
    g = build_duo_grid(loader = (path = "builtin://icosahedral/$level", reader = "builtin_icosahedral"))
    V, F = ESD._subdivide_icosahedron(Float64, level)
    R = g.R
    Nc = ESD.n_cells(g)

    coord = Matrix{Float64}(undef, 3, Nc)   # unit-sphere quadrature point
    area = Vector{Float64}(undef, Nc)
    for c in 1:Nc
        b = corner_binding(V, F[1, c], F[2, c], F[3, c])
        coord[1, c] = ESD.eval_coeff(UCOMP[1], b)
        coord[2, c] = ESD.eval_coeff(UCOMP[2], b)
        coord[3, c] = ESD.eval_coeff(UCOMP[3], b)
        area[c] = g.area[c]
    end

    # Sample the two densities at the quadrature points (single-pathway).
    rho_uniform = ones(Float64, Nc)
    rho_sampled = [ESD.eval_coeff(RHO_AST, Dict{String, Float64}("z" => coord[3, c])) for c in 1:Nc]

    # FAQ integrands for each density.
    function integrands(rho)
        mass = [ESD.eval_coeff(MASS_AST, Dict{String, Float64}("rho" => rho[c], "area" => area[c])) for c in 1:Nc]
        mom = Matrix{Float64}(undef, 3, Nc)
        for c in 1:Nc, d in 1:3
            mom[d, c] = ESD.eval_coeff(MOMENT_AST, Dict{String, Float64}("mass" => mass[c], "coord" => coord[d, c]))
        end
        return mass, mom
    end
    mass_u, mom_u = integrands(rho_uniform)
    mass_s, mom_s = integrands(rho_sampled)

    return Dict{String, Any}(
        "description" => "Canonical cross-binding byte/value contract (bead esd-e5m.1 / D1) for the MPAS-SCVT background quadrature mesh + density sampling FAQ at DUO icosahedral level $level. Every conforming binding (Julia/Rust/Python) MUST reproduce these exact Float64 values, all produced through the single-pathway evaluator: bg_coord (the unit-sphere primal-cell centroid = quadrature point), area (the spherical-triangle weight dA, m^2, from primal_geometry/esd-heg.6), and the density-weighted integrands the SCVT step sums -- bg_mass = rho*area and bg_moment = bg_mass*bg_coord -- for the uniform density (rho==1) and a sampled latitude-graded density (rho = 2 + bg_coord_z). The contract is bitwise Float64 identity (parse then compare bit patterns). Pinned by test/test_mpas_scvt_background_faq.jl. Regenerate with discretizations/grids/mpas/scvt/regenerate_fixtures.jl.",
        "level" => level,
        "n_cells" => Nc,
        "R" => R,
        "density_sample" => "rho(x,y,z) = 2 + z  (z = unit-centroid z-component; rho in [1,3], always positive)",
        "bg_coord_x" => coord[1, :],
        "bg_coord_y" => coord[2, :],
        "bg_coord_z" => coord[3, :],
        "area" => area,
        "rho_uniform" => rho_uniform,
        "bg_mass_uniform" => mass_u,
        "bg_moment_uniform_x" => mom_u[1, :],
        "bg_moment_uniform_y" => mom_u[2, :],
        "bg_moment_uniform_z" => mom_u[3, :],
        "rho_sampled" => rho_sampled,
        "bg_mass_sampled" => mass_s,
        "bg_moment_sampled_x" => mom_s[1, :],
        "bg_moment_sampled_y" => mom_s[2, :],
        "bg_moment_sampled_z" => mom_s[3, :],
    )
end

function main()
    for level in (0, 1)
        d = build_level(level)
        path = joinpath(OUT_DIR, "background_quadrature_level$(level).json")
        write(path, canonical_bytes(d))
        println("wrote $(path) (level=$(level), n_cells=$(d["n_cells"]))")
    end
    return
end

main()
