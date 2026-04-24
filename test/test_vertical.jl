@testsnippet VerticalSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "vertical: sigma uniform from nz" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(coordinate = :sigma, nz = 10)
    @test g isa VerticalGrid
    @test eltype(g) === Float64
    @test ndims(g) == 1
    @test n_cells(g) == 10
    @test n_vertices(g) == 11
    @test n_edges(g) == 10
    # Levels strictly decreasing from 1 → 0.
    @test cell_centers(g, 1) ≈ 0.95
    @test cell_centers(g, 10) ≈ 0.05
    @test all(cell_widths(g) .≈ 0.1)
    @test cell_widths(g, 5) ≈ 0.1
end

@testitem "vertical: sigma from explicit levels" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :sigma, levels = [1.0, 0.7, 0.3, 0.0]
    )
    @test n_cells(g) == 3
    @test cell_widths(g) ≈ [0.3, 0.4, 0.3]
    @test cell_centers(g) ≈ [0.85, 0.5, 0.15]
end

@testitem "vertical: z geometric, strictly increasing" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :z, levels = [0.0, 500.0, 1500.0, 3000.0]
    )
    @test n_cells(g) == 3
    @test cell_widths(g) ≈ [500.0, 1000.0, 1500.0]
    @test metric_eval(g, :z, 2) ≈ 1000.0
    @test metric_eval(g, :dz, 3) ≈ 1500.0
end

@testitem "vertical: theta and z_star require explicit levels" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :theta, nz = 5
    )
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :z_star, nz = 5
    )
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :theta, levels = [280.0, 290.0, 310.0, 350.0]
    )
    @test n_cells(g) == 3
end

@testitem "vertical: eta hybrid via ak/bk" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    # 3-layer eta: bk=1 at surface → 0 aloft, ak ramps up aloft.
    bk = [1.0, 0.7, 0.3, 0.0]
    ak = [0.0, 1.0e4, 3.0e4, 5.0e4]
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :eta, ak = ak, bk = bk, p0 = 1.0e5
    )
    @test n_cells(g) == 3
    # synthesized sigma = ak/p0 + bk should match levels.
    @test g.levels ≈ ak ./ 1.0e5 .+ bk
    # pressure metric averages interface pressures.
    p_lo = ak[1] + bk[1] * 1.0e5
    p_hi = ak[2] + bk[2] * 1.0e5
    @test metric_eval(g, :pressure, 1) ≈ (p_lo + p_hi) / 2
    @test metric_eval(g, :ak, 1) ≈ (ak[1] + ak[2]) / 2
    @test metric_eval(g, :bk, 2) ≈ (bk[2] + bk[3]) / 2
end

@testitem "vertical: eta rejects non-monotone hybrid" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    # Non-decreasing synthesized sigma → DomainError.
    @test_throws DomainError EarthSciDiscretizations.grids.vertical(
        coordinate = :eta,
        ak = [0.0, 0.0, 0.0, 0.0],
        bk = [1.0, 0.5, 0.6, 0.0],
    )
end

@testitem "vertical: hybrid_sigma_theta uniform" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :hybrid_sigma_theta, nz = 8
    )
    @test n_cells(g) == 8
    @test all(cell_widths(g) .≈ 1 / 8)
end

@testitem "vertical: neighbors at ends and middle" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(coordinate = :sigma, nz = 4)
    @test neighbors(g, 1) == Dict(:up => 2)
    @test neighbors(g, 4) == Dict(:down => 3)
    @test neighbors(g, 2) == Dict(:down => 1, :up => 3)
end

@testitem "vertical: bounds and unknown metric" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(coordinate = :sigma, nz = 3)
    @test_throws BoundsError cell_centers(g, 0)
    @test_throws BoundsError cell_centers(g, 4)
    @test_throws BoundsError metric_eval(g, :dz, 0)
    @test_throws ArgumentError metric_eval(g, :nope, 1)
    # :pressure / :ak / :bk on coords without hybrid coefficients.
    @test_throws ArgumentError metric_eval(g, :pressure, 1)
    @test_throws ArgumentError metric_eval(g, :ak, 1)
    @test_throws ArgumentError metric_eval(g, :bk, 1)
end

@testitem "vertical: argument-error contract" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical()
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :unknown_kind, nz = 5
    )
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :sigma, nz = 0
    )
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :sigma, nz = 4, ghosts = -1
    )
    # eta missing ak / bk.
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :eta, bk = [1.0, 0.0]
    )
    @test_throws ArgumentError EarthSciDiscretizations.grids.vertical(
        coordinate = :eta, ak = [0.0, 0.0]
    )
end

@testitem "vertical: to_esm declarative shape" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g = EarthSciDiscretizations.grids.vertical(coordinate = :sigma, nz = 4)
    doc = to_esm(g)
    @test doc["family"] == "vertical"
    @test doc["topology"] == "column"
    @test doc["ndim"] == 1
    @test doc["n_cells"] == 4
    @test doc["n_vertices"] == 5
    @test doc["dtype"] == "float64"
    opts = doc["options"]
    @test opts["coordinate"] == "sigma"
    @test opts["nz"] == 4
    @test length(opts["levels"]) == 5
    # No hybrid coeffs for plain sigma → keys absent.
    @test !haskey(opts, "ak")
    @test !haskey(opts, "bk")
    prov = doc["provenance"]
    @test prov["family"] == "vertical"
    @test prov["binding"] == "julia"
    @test prov["coordinate"] == "sigma"
end

@testitem "vertical: to_esm carries hybrid coefficients for eta" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    bk = [1.0, 0.4, 0.0]
    ak = [0.0, 2.0e4, 5.0e4]
    g = EarthSciDiscretizations.grids.vertical(
        coordinate = :eta, ak = ak, bk = bk, p0 = 1.0e5
    )
    doc = to_esm(g)
    @test doc["options"]["ak"] == ak
    @test doc["options"]["bk"] == bk
    @test doc["options"]["p0"] ≈ 1.0e5
end

@testitem "vertical: family and dtype propagation" setup = [VerticalSetup] tags = [:grid, :vertical] begin
    g32 = EarthSciDiscretizations.grids.vertical(
        coordinate = :sigma, nz = 4, dtype = Float32
    )
    @test eltype(g32) === Float32
    @test g32.dtype == "float32"
    @test family(g32) == "vertical"
end
