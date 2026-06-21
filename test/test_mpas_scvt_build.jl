# MPAS-SCVT mesh generator host DRIVER (esd-e5m.4 / D4): the external Lloyd
# fixed-point loop over the declarative D2 step (lloyd_step.esm, via the
# value-invention front-door materialize_value_invention), capped ONCE at
# convergence by the D3 spherical-topology leaf (scvt_voronoi_connectivity),
# emitting the same MpasMeshData the imperative DUO Voronoi dual produces.
#
# Per the epic's DECLARATIVE-OR-FAIL contract the per-iteration STEP is declarative
# and the LOOP is host RHS-only (the convergence test + the sphere re-projection,
# which needs a sqrt outside the value-invention op set). These tests prove the
# driver: (1) the loop converges, (2) it emits a valid mass-conserving MpasMeshData
# with the icosahedral-dual topology, (3) a density function drives variable
# resolution, (4) the output is a genuine fixed point of ONE more declarative step
# (the host loop wraps the declarative step — no recurrence in the IR).

@testitem "scvt driver: external loop → dodecahedral MpasMeshData (uniform CVT)" tags = [:grid, :mpas, :scvt, :driver, :lloyd] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using LinearAlgebra: norm

    R = 6.371e6
    sphere_area = 4 * pi * R^2

    # Seed = the 12 level-0 icosahedral vertices; background FINER than the seed so
    # every generator is attended (LLOYD_STEP_CONTRACT.md §4).
    V0, _ = ESD.duo_subdivide_faq(Float64, 0)               # (3, 12)
    @test size(V0) == (3, 12)

    # The host fixed-point loop converges (the icosahedral generators are a CVT).
    bg_coord, bg_mass = scvt_background_quadrature(2)
    @test size(bg_coord) == (320, 3) && length(bg_mass) == 320
    @test sum(bg_mass) ≈ sphere_area rtol = 1e-12            # D1 mass = ∫ρ dA = 4πR²
    sol = scvt_lloyd_solve(V0, bg_coord, bg_mass; tol = 1e-10, max_iters = 200)
    @test sol.converged
    @test sol.iterations >= 1
    @test size(sol.generators) == (3, 12)
    @test all(c -> isapprox(norm(sol.generators[:, c]), R; rtol = 1e-12), 1:12)

    # The full driver emits the dodecahedron (dual of the icosahedron): 12 pentagons.
    mesh = build_scvt_mesh(; generators = V0, density = nothing, background_level = 2,
        tol = 1e-10, max_iters = 200)
    @test mesh isa MpasMeshData
    @test mesh.n_cells == 12
    @test mesh.n_edges == 30                                  # 3n-6
    @test mesh.n_vertices == 20                               # 2n-4 Delaunay triangles
    @test mesh.max_edges == 5
    @test all(==(5), mesh.n_edges_on_cell)                    # every cell a pentagon

    # Mass conservation: the Voronoi cells tile the sphere.
    @test sum(mesh.area_cell) ≈ sphere_area rtol = 1e-12
    @test all(>(0), mesh.area_cell)

    # Cell centres lie on the sphere of radius R and lon/lat agree with the
    # cartesian. The lon/lat↔cartesian round-trip is compared with an R-scaled
    # ABSOLUTE tolerance (a relative tol trips on a near-zero coordinate component,
    # e.g. a generator on the y-z plane has x_cell ≈ 0).
    for c in 1:12
        @test norm((mesh.x_cell[c], mesh.y_cell[c], mesh.z_cell[c])) ≈ R rtol = 1e-12
        @test mesh.z_cell[c] ≈ R * sin(mesh.lat_cell[c]) atol = 1e-6 * R
        @test mesh.x_cell[c] ≈ R * cos(mesh.lat_cell[c]) * cos(mesh.lon_cell[c]) atol = 1e-6 * R
        @test mesh.y_cell[c] ≈ R * cos(mesh.lat_cell[c]) * sin(mesh.lon_cell[c]) atol = 1e-6 * R
    end

    # Edge arrays are consistent: each edge joins its two cells, both arcs positive.
    @test size(mesh.cells_on_edge) == (2, 30)
    @test all(1 .<= mesh.cells_on_edge .<= 12)
    @test all(>(0), mesh.dc_edge) && all(>(0), mesh.dv_edge)

    # The adjacency is a valid closed mesh (bounds + neighbour reciprocity).
    @test ESD.check_mesh(mesh, true) === nothing
end

@testitem "scvt driver: icosahedral-dual at level 1 (pentagons + hexagons)" tags = [:grid, :mpas, :scvt, :driver] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations

    R = 6.371e6
    V1, _ = ESD.duo_subdivide_faq(Float64, 1)               # (3, 42)
    mesh = build_scvt_mesh(; generators = V1, density = nothing, background_level = 3,
        tol = 1e-9, max_iters = 300)
    @test mesh.n_cells == 42
    @test mesh.n_edges == 120                                # 3·42-6
    @test mesh.n_vertices == 80                              # 2·42-4
    @test mesh.max_edges == 6
    # The icosahedral dual has exactly 12 pentagons; the rest are hexagons.
    @test count(==(5), mesh.n_edges_on_cell) == 12
    @test count(==(6), mesh.n_edges_on_cell) == 30
    @test sum(mesh.area_cell) ≈ 4 * pi * R^2 rtol = 1e-12
    @test ESD.check_mesh(mesh, true) === nothing
end

@testitem "scvt driver: a density function drives variable resolution" tags = [:grid, :mpas, :scvt, :driver, :lloyd] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations

    R = 6.371e6
    V0, _ = ESD.duo_subdivide_faq(Float64, 0)

    mesh_u = build_scvt_mesh(; generators = V0, density = nothing, background_level = 2)
    # ρ = 2 + z (z = unit-centroid z; ρ ∈ [1, 3], always positive) — the same
    # smooth latitude-graded density D1/D2 pin. It reweights the centroids, so the
    # converged mesh (and its cell areas) differ from the uniform CVT.
    mesh_d = build_scvt_mesh(; generators = V0, density = (x, y, z) -> 2.0 + z,
        background_level = 2)

    @test mesh_d.area_cell != mesh_u.area_cell               # density changed the result
    @test mesh_d.lat_cell != mesh_u.lat_cell
    # Still a valid, mass-conserving tiling of the whole sphere.
    @test mesh_d.n_cells == 12
    @test sum(mesh_d.area_cell) ≈ 4 * pi * R^2 rtol = 1e-12
    @test ESD.check_mesh(mesh_d, true) === nothing
end

@testitem "scvt driver: the host loop wraps the declarative step (genuine fixed point)" tags = [:grid, :mpas, :scvt, :driver, :lloyd] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using LinearAlgebra: norm

    # A deterministic OFF-CVT seed: nudge each icosahedral vertex by a fixed pattern
    # then renormalise. The Lloyd loop must do real work (> 1 iteration) and pull it
    # back to a centroidal Voronoi fixed point.
    V0, _ = ESD.duo_subdivide_faq(Float64, 0)
    Vp = copy(V0)
    for i in 1:size(Vp, 2)
        Vp[:, i] .+= 0.15 .* Float64[(i % 3) - 1, ((i + 1) % 3) - 1, ((i + 2) % 3) - 1]
        Vp[:, i] ./= norm(Vp[:, i])
    end

    bg = scvt_background_quadrature(2)
    sol = scvt_lloyd_solve(Vp, bg...; tol = 1e-10, max_iters = 500)
    @test sol.converged
    @test sol.iterations > 1                                  # the seed was off the fixed point

    # The converged generators are a FIXED POINT of ONE more declarative step: the
    # host loop iterates the SAME declarative Lloyd step, so applying it once more
    # barely moves them. This is the behavioural proof that the loop is host
    # RHS-only over the declarative step — no recurrence is lowered into the IR.
    sol_unit = Float64[sol.generators[d, c] / 6.371e6 for d in 1:3, c in 1:size(sol.generators, 2)]
    step1 = scvt_lloyd_solve(sol_unit, bg...; tol = 1e-14, max_iters = 1)
    @test step1.displacement < 1e-9

    # And the driver from the same perturbed seed still emits a valid mesh.
    mesh = build_scvt_mesh(; generators = Vp, density = nothing, background_level = 2,
        tol = 1e-10, max_iters = 500)
    @test ESD.check_mesh(mesh, true) === nothing
    @test sum(mesh.area_cell) ≈ 4 * pi * 6.371e6^2 rtol = 1e-12
end

@testitem "scvt driver: input validation" tags = [:grid, :mpas, :scvt, :driver] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations

    V0, _ = ESD.duo_subdivide_faq(Float64, 0)
    # Fewer than 4 generators cannot form a closed spherical mesh.
    @test_throws ArgumentError build_scvt_mesh(; generators = V0[:, 1:3], background_level = 1)
    # Generators must be (3, Nc) Cartesian columns.
    @test_throws ArgumentError build_scvt_mesh(; generators = permutedims(V0), background_level = 1)
    @test_throws ArgumentError scvt_lloyd_solve(permutedims(V0), scvt_background_quadrature(1)...)
    # Degenerate parameters are rejected.
    @test_throws DomainError scvt_background_quadrature(-1)
    @test_throws DomainError scvt_lloyd_solve(V0, scvt_background_quadrature(1)...; tol = 0.0)
end
