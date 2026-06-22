# MPAS-SCVT conformance DRAIN (esd-e5m.5 / D5) — the driver-level conformance of
# the declarative MPAS-SCVT mesh generator (epic esd-e5m). The lower layers are
# already pinned: D1 background quadrature (test_mpas_scvt_background_faq.jl), D2
# Lloyd step (test_mpas_scvt_lloyd_step_faq.jl), D3 spherical-topology leaf +
# cross-binding determinism golden (test_mpas_scvt_topology_leaf.jl). This file
# drains the remaining DRIVER-level conformance obligations:
#
#   (a) ρ≡1 regression — the full Lloyd driver (build_scvt_mesh) recovers the
#       quasi-uniform icosahedral DUO-dual ("x1.*" family, _duo_voronoi_dual): the
#       same mesh structure, and a geometry that APPROACHES the DUO-dual as the
#       background quadrature refines (the rigorous "reproduces to tolerance").
#   (b) CVT fixed-point property — at convergence each generator is the
#       density-weighted centroid of its Voronoi cell (uniform ρ≡1 AND a variable
#       density ρ=2+z), checked both through the declarative D2 step and an
#       independent test-local recompute.
#   (c) variable-resolution reference fixture — the converged ρ=2+z mesh matches
#       the checked-in internal reference golden to tolerance, with the
#       density-shift correctness property (smaller cells where ρ is larger).
#
# Cross-binding: the DETERMINISTIC integer topology is binding-neutral and pinned
# by the D3 leaf golden (tests/conformance/grids/mpas/scvt/topology_leaf/); the
# TOLERANCE geometry rides the shared duo_dual_geometry_faq §5.8 contract. The
# Lloyd LOOP itself is the host driver (reference binding Julia) — see
# discretizations/grids/mpas/scvt/BUILD_DRIVER_CONTRACT.md and
# tests/conformance/grids/mpas/scvt/variable_resolution/README.md.

@testitem "scvt conformance (a): ρ≡1 driver recovers the quasi-uniform DUO-dual under refinement" tags = [:conformance, :grid, :mpas, :scvt, :driver, :regression] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations

    R = 6.371e6
    V1, _ = ESD.duo_subdivide_faq(Float64, 1)              # 42-generator level-1 seed
    imp = ESD._duo_voronoi_dual(1; R = R)                  # the DUO-dual "x1.*" reference

    # Full Lloyd driver at two background refinements (uniform density → CVT).
    m_coarse = build_scvt_mesh(;
        generators = V1, density = nothing,
        background_level = 2, R = R, tol = 1.0e-10, max_iters = 2000
    )
    m_fine = build_scvt_mesh(;
        generators = V1, density = nothing,
        background_level = 3, R = R, tol = 1.0e-10, max_iters = 2000
    )

    # (1) SAME mesh family as the DUO-dual — exact integer structure at every level.
    for m in (m_coarse, m_fine)
        @test m.n_cells == imp.n_cells == 42
        @test m.n_edges == imp.n_edges == 120              # 3n−6
        @test m.n_vertices == imp.n_vertices == 80         # 2n−4
        @test m.max_edges == imp.max_edges == 6
        # Identical valence histogram: 12 pentagons + 30 hexagons (icosahedral dual).
        @test count(==(5), m.n_edges_on_cell) == 12
        @test count(==(6), m.n_edges_on_cell) == 30
        @test sort(m.n_edges_on_cell) == sort(imp.n_edges_on_cell)
        # Mass conservation: the Voronoi cells tile the whole sphere.
        @test sum(m.area_cell) ≈ 4 * pi * R^2 rtol = 1.0e-12
        @test all(>(0), m.area_cell)
        @test ESD.check_mesh(m, true) === nothing
    end

    # (2) Quasi-uniform: the icosahedral-dual family keeps cell areas within a tight
    # band. The DISCRETE CVT over a finite background is slightly more dispersed than
    # the smooth raw geodesic DUO-dual (its area max/min ≈ 1.22 at bg level 2 vs the
    # DUO-dual's ≈ 1.13) and tightens toward it as the background refines — that
    # convergence is carried by the area gap in (3); here we just pin that both stay
    # quasi-uniform.
    @test maximum(imp.area_cell) / minimum(imp.area_cell) < 1.3
    for m in (m_coarse, m_fine)
        @test maximum(m.area_cell) / minimum(m.area_cell) < 1.3
    end

    # (3) Refinement convergence — the headline ρ≡1 regression. The driver geometry
    # APPROACHES the DUO-dual as the background quadrature refines: the per-cell area
    # gap shrinks monotonically (the discrete-CVT discretization error → 0). This is
    # the rigorous form of "uniform-density SCVT reproduces the x1.* family to
    # tolerance" — the tolerance is the background quadrature resolution.
    gap(m) = maximum(abs.(m.area_cell .- imp.area_cell) ./ imp.area_cell)
    @test gap(m_fine) < gap(m_coarse)                      # refinement closes the gap
    @test gap(m_fine) < 0.08                               # within tolerance at bg level 3

    # The CONNECTIVITY recovers the DUO-dual's exactly when the generators are NOT
    # Lloyd-moved (fed straight to the leaf) — the byte-identical ρ≡1 topology
    # regression pinned in test_mpas_scvt_topology_leaf.jl. The full driver moves the
    # generators to the centroidal fixed point, so a few near-cocircular edges may
    # flip; the VALENCE histogram (asserted above) is the topology invariant that
    # survives the move.
    leaf = scvt_voronoi_connectivity(V1; R = R)
    @test leaf.cells_on_cell == imp.cells_on_cell          # leaf-on-seed = DUO-dual
end

@testitem "scvt conformance (b): CVT property — each generator is its cell's density-weighted centroid" tags = [:conformance, :grid, :mpas, :scvt, :driver, :lloyd, :cvt] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using LinearAlgebra: norm

    R = 6.371e6
    V1, _ = ESD.duo_subdivide_faq(Float64, 1)

    # An independent, test-local Voronoi-centroid oracle (the conformance pattern:
    # re-derive the expected value without the driver's machinery). For converged
    # generators `gu` (unit) over the background (`bgc` unit points, `bgm` density-
    # weighted measure), assign each point to its nearest generator (squared
    # Euclidean, smallest-id tie-break — the §5.7 rule-6 convention) and form the
    # density-weighted centroid; the CVT property is `normalize(centroid_g) == gu_g`.
    function centroid_residual(gu, bgc, bgm)
        Nc = size(gu, 2)
        num = zeros(3, Nc)
        den = zeros(Nc)
        for p in 1:size(bgc, 1)
            best = 0
            bd = Inf
            for g in 1:Nc
                d = (bgc[p, 1] - gu[1, g])^2 + (bgc[p, 2] - gu[2, g])^2 + (bgc[p, 3] - gu[3, g])^2
                if d < bd            # strict `<` keeps the smallest-id generator on ties
                    bd = d
                    best = g
                end
            end
            num[1, best] += bgm[p] * bgc[p, 1]
            num[2, best] += bgm[p] * bgc[p, 2]
            num[3, best] += bgm[p] * bgc[p, 3]
            den[best] += bgm[p]
        end
        maxr = 0.0
        n_attended = 0
        for g in 1:Nc
            den[g] == 0 && continue          # an unattended generator is held, not centroidal
            n_attended += 1
            c = num[:, g] ./ den[g]
            c ./= norm(c)
            maxr = max(maxr, norm(c .- gu[:, g]))
        end
        return maxr, n_attended
    end

    for density in (nothing, (x, y, z) -> 2.0 + z)         # uniform CVT and variable-resolution SCVT
        bgc, bgm = scvt_background_quadrature(3; density = density, R = R)
        sol = scvt_lloyd_solve(V1, bgc, bgm; R = R, tol = 1.0e-12, max_iters = 2000)
        @test sol.converged
        gu = sol.generators ./ R                            # converged generators, unit sphere

        # Through the DECLARATIVE D2 step: one more Lloyd iteration on the converged
        # generators is `max_g ‖normalize(centroid_g) − gen_g‖` — the CVT residual.
        # ≈ 0 means each generator already IS its cell's density-weighted centroid.
        step1 = scvt_lloyd_solve(gu, bgc, bgm; R = R, tol = 1.0e-16, max_iters = 1)
        @test step1.displacement < 1.0e-9

        # And against the independent test-local centroid oracle (summation order
        # differs from the relational group-aggregate, hence a looser threshold). The
        # background (level 3) is finer than the 42 generators, so every generator is
        # attended (none held at its seed) — the CVT property covers the whole mesh.
        resid, n_attended = centroid_residual(gu, bgc, bgm)
        @test n_attended == 42
        @test resid < 1.0e-6
    end
end

@testitem "scvt conformance (c): variable-resolution reference fixture (ρ=2+z)" tags = [:conformance, :grid, :mpas, :scvt, :driver, :varres] begin
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON

    # Pearson correlation, inlined (Statistics is not a dependency of this package).
    function _corr(a, b)
        ma = sum(a) / length(a)
        mb = sum(b) / length(b)
        da = a .- ma
        db = b .- mb
        return sum(da .* db) / sqrt(sum(da .^ 2) * sum(db .^ 2))
    end

    R = 6.371e6
    GOLDEN = JSON.parsefile(
        joinpath(
            @__DIR__, "..", "tests", "conformance", "grids",
            "mpas", "scvt", "variable_resolution", "golden.json"
        )
    )
    gref = GOLDEN["variable_resolution"]
    uref = GOLDEN["uniform_density"]

    V1, _ = ESD.duo_subdivide_faq(Float64, 1)
    mesh = build_scvt_mesh(;
        generators = V1, density = (x, y, z) -> 2.0 + z,
        background_level = 3, R = R, tol = 1.0e-12, max_iters = 2000
    )

    # Exact integer structure.
    @test mesh.n_cells == gref["n_cells"] == 42
    @test mesh.n_edges == gref["n_edges"] == 120
    @test mesh.n_vertices == gref["n_vertices"] == 80
    @test mesh.max_edges == gref["max_edges"] == 6
    @test mesh.n_edges_on_cell == Int.(gref["n_edges_on_cell"])

    # Geometry to the §5.8 tolerance contract (the reference is bitwise reproducible,
    # so the tolerance is generous margin, not slop).
    approx_vec(a, b; rtol = 1.0e-8, atol = 1.0e-6) = all(isapprox.(a, Float64.(b); rtol = rtol, atol = atol))
    @test approx_vec(mesh.area_cell, gref["area_cell"]; atol = 1.0e-3)   # m² — atol vs R²-scale
    @test approx_vec(mesh.lat_cell, gref["lat_cell"])
    @test approx_vec(mesh.lon_cell, gref["lon_cell"])
    @test approx_vec(mesh.x_cell, gref["x_cell"]; atol = 1.0e-3)
    @test approx_vec(mesh.y_cell, gref["y_cell"]; atol = 1.0e-3)
    @test approx_vec(mesh.z_cell, gref["z_cell"]; atol = 1.0e-3)

    # Correctness — the variable-resolution SIGNAL: ρ = 2 + z is largest toward the
    # north pole (z → 1), so the converged mesh has SMALLER cells there — area is
    # negatively correlated with z.
    @test _corr(mesh.area_cell, mesh.z_cell) < -0.3

    # The density genuinely reshaped the mesh: the var-res areas differ materially
    # from the uniform-density (ρ≡1 CVT) companion pinned alongside (≈ 17% at the
    # most-shifted cell), and the live var-res rebuild differs from that companion.
    @test Float64.(gref["area_cell"]) != Float64.(uref["area_cell"])
    @test maximum(abs.(mesh.area_cell .- Float64.(uref["area_cell"])) ./ Float64.(uref["area_cell"])) > 0.05

    # Still a valid, mass-conserving tiling of the whole sphere.
    @test sum(mesh.area_cell) ≈ 4 * pi * R^2 rtol = 1.0e-12
    @test all(>(0), mesh.area_cell)
    @test ESD.check_mesh(mesh, true) === nothing
end
