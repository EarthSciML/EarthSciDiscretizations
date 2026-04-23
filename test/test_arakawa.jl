@testsnippet ArakawaSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Arakawa stagger -> variable locations" setup = [ArakawaSetup] tags = [:arakawa] begin
    # A: all colocated at cell centers.
    @test arakawa_variable_locations(ArakawaA) === (CellCenter, CellCenter, CellCenter)
    # B: h at center, u/v at corners.
    @test arakawa_variable_locations(ArakawaB) === (CellCenter, Corner, Corner)
    # C: h center, u u-face, v v-face.
    @test arakawa_variable_locations(ArakawaC) === (CellCenter, UEdge, VEdge)
    # D: h center, u v-face, v u-face (swapped).
    @test arakawa_variable_locations(ArakawaD) === (CellCenter, VEdge, UEdge)
    # E: topologically like B (rotated 45°; rotation flag carried separately).
    @test arakawa_variable_locations(ArakawaE) === (CellCenter, Corner, Corner)
end

@testitem "Arakawa location shape" setup = [ArakawaSetup] tags = [:arakawa] begin
    @test arakawa_location_shape(CellCenter, 10, 20) == (10, 20)
    @test arakawa_location_shape(UEdge, 10, 20) == (11, 20)
    @test arakawa_location_shape(VEdge, 10, 20) == (10, 21)
    @test arakawa_location_shape(Corner, 10, 20) == (11, 21)
end

@testitem "CartesianBase constructor validates inputs" setup = [ArakawaSetup] tags = [:arakawa] begin
    @test_throws DomainError CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 0, ny = 4)
    @test_throws DomainError CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = -1)
    @test_throws DomainError CartesianBase(xlo = 1.0, xhi = 0.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    @test_throws DomainError CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 1.0, yhi = 1.0, nx = 4, ny = 4)
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    @test b.nx == 4 && b.ny == 4
end

@testitem "ArakawaGrid keyword generator" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 2.0, nx = 10, ny = 20)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)
    @test g isa ArakawaGrid
    @test g.stagger === ArakawaC
    @test g.ghosts == 0
    @test eltype(g) === Float64

    # ArakawaStagger value form works too.
    g2 = EarthSciDiscretizations.grids.arakawa(base = b, stagger = ArakawaB, ghosts = 3, dtype = Float32)
    @test g2.stagger === ArakawaB
    @test g2.ghosts == 3
    @test eltype(g2) === Float32
end

@testitem "ArakawaGrid generator error contract" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    # Missing required options → ArgumentError (per GRIDS_API.md §9, Julia row).
    @test_throws ArgumentError EarthSciDiscretizations.grids.arakawa(stagger = :C)
    @test_throws ArgumentError EarthSciDiscretizations.grids.arakawa(base = b)
    # Invalid option value → DomainError.
    @test_throws DomainError EarthSciDiscretizations.grids.arakawa(base = b, stagger = :Q)
end

@testitem "Arakawa C-grid accessors recover Cartesian coordinates" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 10, ny = 10)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)

    # Cell center (1,1) sits at half-cell from the lower corner.
    cx, cy = cell_centers(g, 1, 1)
    @test cx ≈ 0.05
    @test cy ≈ 0.05

    # Cell center (10,10) is at the far upper-right half-cell.
    cx, cy = cell_centers(g, 10, 10)
    @test cx ≈ 0.95
    @test cy ≈ 0.95

    # For C: u lives at u-faces. u(1,1) is on the western domain boundary at
    # mid-y of row 1.
    ux, uy = u_face(g, 1, 1)
    @test ux ≈ 0.0
    @test uy ≈ 0.05
    # u(11, 1) is the eastern boundary of the last column.
    ux, uy = u_face(g, 11, 1)
    @test ux ≈ 1.0
    @test uy ≈ 0.05

    # For C: v lives at v-faces. v(1,1) is on the southern boundary at mid-x.
    vx, vy = v_face(g, 1, 1)
    @test vx ≈ 0.05
    @test vy ≈ 0.0
    vx, vy = v_face(g, 1, 11)
    @test vx ≈ 0.05
    @test vy ≈ 1.0

    # Corners
    @test corners(g, 1, 1) == (0.0, 0.0)
    @test corners(g, 11, 11) == (1.0, 1.0)
end

@testitem "Arakawa A-grid colocates u and v at cell centers" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :A)
    # h, u, v all live at (nx, ny) cell-centers.
    @test variable_shape(g, :h) == (4, 4)
    @test variable_shape(g, :u) == (4, 4)
    @test variable_shape(g, :v) == (4, 4)
    # u and v at (2,3) equal the cell center there.
    @test u_face(g, 2, 3) == cell_centers(g, 2, 3)
    @test v_face(g, 2, 3) == cell_centers(g, 2, 3)
end

@testitem "Arakawa B-grid puts u,v at corners" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :B)
    @test variable_shape(g, :h) == (4, 4)
    @test variable_shape(g, :u) == (5, 5)
    @test variable_shape(g, :v) == (5, 5)
    @test u_face(g, 1, 1) == corners(g, 1, 1)
    @test v_face(g, 3, 2) == corners(g, 3, 2)
end

@testitem "Arakawa D-grid swaps C's u/v faces" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 2.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    gC = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)
    gD = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :D)
    # D's u lives where C's v lives; D's v where C's u lives.
    @test variable_shape(gD, :u) == variable_shape(gC, :v)
    @test variable_shape(gD, :v) == variable_shape(gC, :u)
    @test u_face(gD, 1, 1) == v_face(gC, 1, 1)
    @test v_face(gD, 1, 1) == u_face(gC, 1, 1)
end

@testitem "Arakawa neighbors bound-check and boundary handling" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)

    # Interior cell has all four neighbors.
    w, e, s, n = neighbors(g, 2, 2)
    @test w === (CellCenter, 1, 2)
    @test e === (CellCenter, 3, 2)
    @test s === (CellCenter, 2, 1)
    @test n === (CellCenter, 2, 3)

    # Corner cells lose two sides.
    w, e, s, n = neighbors(g, 1, 1)
    @test w === nothing
    @test s === nothing
    @test e === (CellCenter, 2, 1)
    @test n === (CellCenter, 1, 2)

    # Bounds check.
    @test_throws DomainError neighbors(g, 0, 1)
    @test_throws DomainError neighbors(g, 1, 5)

    # UEdge neighbors use the UEdge index range (nx+1, ny).
    w, e, s, n = neighbors(g, UEdge, 5, 4)  # last column, last row
    @test e === nothing
    @test n === nothing
    @test w === (UEdge, 4, 4)
    @test s === (UEdge, 5, 3)
end

@testitem "Arakawa metric_eval yields consistent Cartesian spacing" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 2.0, ylo = 0.0, yhi = 6.0, nx = 4, ny = 3)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C)
    # dx = 2/4 = 0.5; dy = 6/3 = 2.0; area = 1.0.
    @test metric_eval(g, :dx, 1, 1) ≈ 0.5
    @test metric_eval(g, :dy, 1, 1) ≈ 2.0
    @test metric_eval(g, :area, 4, 3) ≈ 1.0
    @test_throws ArgumentError metric_eval(g, :not_a_metric, 1, 1)
    @test_throws DomainError metric_eval(g, :dx, 5, 1)
end

@testitem "Arakawa dtype propagates through metrics" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    g32 = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C, dtype = Float32)
    @test metric_eval(g32, :dx, 1, 1) isa Float32
end

@testitem "Arakawa to_esm produces schema-valid declarative config" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 5, ny = 5)
    g = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :C, ghosts = 2)
    d = to_esm(g)

    @test d isa Dict{String, Any}
    @test d["family"] == "arakawa"
    @test d["stagger"] == "C"
    @test d["dtype"] == "float64"
    @test d["topology"] == "block_structured"
    @test d["ghosts"] == 2
    @test d["n_cells"] == 25
    @test d["rotated"] == false

    # Per §0 correction: no inline geometry arrays. The base is a declarative ref.
    @test !haskey(d, "cells")
    @test !haskey(d, "edges")
    @test !haskey(d, "vertices")

    base = d["base"]
    @test base["family"] == "cartesian"
    @test base["nx"] == 5
    @test base["ny"] == 5
    @test base["extent"] == [[0.0, 0.0], [1.0, 1.0]]
end

@testitem "Arakawa E carries a rotated flag" setup = [ArakawaSetup] tags = [:arakawa] begin
    b = CartesianBase(xlo = 0.0, xhi = 1.0, ylo = 0.0, yhi = 1.0, nx = 4, ny = 4)
    gE = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :E)
    gB = EarthSciDiscretizations.grids.arakawa(base = b, stagger = :B)
    dE = to_esm(gE)
    dB = to_esm(gB)
    @test dE["stagger"] == "E"
    @test dE["rotated"] == true
    @test dB["rotated"] == false
    # Topologically E mirrors B (same location shapes per variable).
    @test variable_shape(gE, :h) == variable_shape(gB, :h)
    @test variable_shape(gE, :u) == variable_shape(gB, :u)
    @test variable_shape(gE, :v) == variable_shape(gB, :v)
end
