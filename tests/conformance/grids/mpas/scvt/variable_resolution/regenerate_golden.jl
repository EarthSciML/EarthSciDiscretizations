#!/usr/bin/env julia
# Regenerate the MPAS-SCVT VARIABLE-RESOLUTION reference fixture
# (epic esd-e5m, bead esd-e5m.5 / D5, conformance criterion (c)).
#
# The headline capability of the declarative MPAS-SCVT generator over the fixed
# icosahedral DUO dual is VARIABLE RESOLUTION via a density function ρ: the Lloyd
# loop pulls more (smaller) cells into the region where ρ is larger. Criterion (c)
# of the D5 drain calls for a checked-in variable-resolution reference fixture.
#
# The standard external references (jigsaw / MPAS-Tools `MpasMeshConverter`) are
# offline mesh generators NOT available in this environment, so — exactly as the
# epic's DECLARATIVE-OR-FAIL culture treats an unavailable external tool — this
# fixture is an INTERNAL reference: the converged variable-resolution mesh that
# the canonical `build_scvt_mesh` pipeline produces, pinned so the var-res result
# is reproducible and its density-shift behaviour is locked against regression.
# The mesh is generated THROUGH THE SINGLE PATHWAY (the D1 background-quadrature
# FAQ + the D2 declarative Lloyd step via `materialize_value_invention`, driven by
# the D4 host loop), per AGENTS.md ("golden regeneration drives the canonical
# pipeline") — no shadow generator.
#
# Canonical instance: the 42-generator level-1 icosahedral seed under the smooth
# latitude-graded density ρ = 2 + z (z = unit-sphere z-component; ρ ∈ [1, 3], the
# same density D1/D2/D4 pin), background level 3 (1280 quadrature points, finer
# than the 42 generators so every generator is attended), converged to a discrete
# centroidal fixed point (tol 1e-12). ρ is larger toward the north pole (z → 1),
# so the converged mesh has SMALLER cells there — the variable-resolution signal
# this fixture pins.
#
# The mesh GEOMETRY (cell areas, geographic coordinates, generator positions)
# rides the §5.8 TOLERANCE contract (it derives from the libm-dependent
# spherical-excess areas and the centroid sums), not byte identity — so the test
# compares to tolerance. The Julia reference is itself bitwise reproducible run to
# run (the discrete Lloyd fixed point is deterministic), so this golden is stable.
#
# Usage (from the repo root):
#     julia --project=. tests/conformance/grids/mpas/scvt/variable_resolution/regenerate_golden.jl

using EarthSciDiscretizations
const ESD = EarthSciDiscretizations
using JSON

const R_EARTH = 6.371e6
const SEED_LEVEL = 1            # 42-generator icosahedral seed
const BG_LEVEL = 3             # 1280 background quadrature points (attends all 42)
const TOL = 1e-12
const MAX_ITERS = 2000

# The canonical variable-resolution density: smooth, latitude-graded, ρ ∈ [1, 3].
density(x, y, z) = 2.0 + z

function build_reference()
    V, _ = ESD.duo_subdivide_faq(Float64, SEED_LEVEL)        # (3, 42) seed generators
    mesh = build_scvt_mesh(; generators = V, density = density,
        background_level = BG_LEVEL, R = R_EARTH, tol = TOL, max_iters = MAX_ITERS)
    # Uniform-density (ρ ≡ 1 → CVT) companion: the same seed/background WITHOUT the
    # density, pinned alongside so the test can assert the density genuinely moved
    # the mesh (var-res ≠ uniform) against fixed reference numbers, not just "≠".
    mesh_u = build_scvt_mesh(; generators = V, density = nothing,
        background_level = BG_LEVEL, R = R_EARTH, tol = TOL, max_iters = MAX_ITERS)
    return mesh, mesh_u
end

function mesh_dict(mesh)
    return Dict(
        "n_cells" => mesh.n_cells,
        "n_edges" => mesh.n_edges,
        "n_vertices" => mesh.n_vertices,
        "max_edges" => mesh.max_edges,
        "n_edges_on_cell" => mesh.n_edges_on_cell,
        "area_cell" => mesh.area_cell,
        "lat_cell" => mesh.lat_cell,
        "lon_cell" => mesh.lon_cell,
        "x_cell" => mesh.x_cell,
        "y_cell" => mesh.y_cell,
        "z_cell" => mesh.z_cell,
    )
end

function main()
    mesh, mesh_u = build_reference()
    out = Dict(
        "description" =>
            "MPAS-SCVT VARIABLE-RESOLUTION reference fixture (epic esd-e5m, bead esd-e5m.5 / D5, " *
            "criterion (c)). The converged variable-resolution mesh `build_scvt_mesh` produces from " *
            "the 42-generator level-1 icosahedral seed under the latitude-graded density ρ = 2 + z " *
            "(z = unit-sphere z; ρ ∈ [1,3]), background level 3, tol 1e-12 — generated through the " *
            "canonical single pathway (D1 background-quadrature FAQ + D2 declarative Lloyd step via " *
            "materialize_value_invention, D4 host loop). The mesh GEOMETRY rides the §5.8 tolerance " *
            "contract (compared to tolerance, not byte-identity); the integer structure is exact. " *
            "An INTERNAL reference: the external jigsaw / MPAS-Tools generators are offline tools not " *
            "available in this environment (see VERDICT in README.md), so this pins the canonical " *
            "pipeline's own reproducible var-res output and locks the density-shift behaviour " *
            "(smaller cells toward the dense north pole). The uniform-density (ρ≡1 CVT) companion is " *
            "pinned alongside so the test asserts the density genuinely reshaped the mesh.",
        "reference_binding" => "julia",
        "seed_level" => SEED_LEVEL,
        "background_level" => BG_LEVEL,
        "R" => R_EARTH,
        "tol" => TOL,
        "density" => "rho(x,y,z) = 2 + z",
        "variable_resolution" => mesh_dict(mesh),
        "uniform_density" => mesh_dict(mesh_u),
    )
    dir = @__DIR__
    path = joinpath(dir, "golden.json")
    open(path, "w") do io
        JSON.print(io, out, 2)
        write(io, "\n")
    end
    println("wrote ", path)
end

main()
