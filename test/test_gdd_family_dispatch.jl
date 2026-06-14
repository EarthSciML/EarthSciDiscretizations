using Test
using TestItems

# Unit tests for the grid-family dispatch router in _inject_grids! / _inject_rules!.
# Tests that each family selects the right ESM grid block and lowering path.

@testitem "GDD family dispatch: default family is cartesian when omitted" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            # no "family" field → must default to "cartesian"
            "spatial" => Dict{String,Any}(
                "x" => Dict{String,Any}("min" => 0.0, "max" => 1.0, "grid_spacing" => 0.25),
            ),
            "boundary_conditions" => Any[Dict{String,Any}("type" => "periodic", "dimensions" => ["x"])],
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "cartesian"
    @test length(block["dimensions"]) == 1
    @test block["dimensions"][1]["name"] == "x"
    @test block["dimensions"][1]["size"] == 4
    @test block["dimensions"][1]["periodic"] == true
end

@testitem "GDD family dispatch: explicit cartesian family passes through" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "cartesian",
            "spatial" => Dict{String,Any}(
                "x" => Dict{String,Any}("min" => 0.0, "max" => 1.0, "grid_spacing" => 0.1),
            ),
            "boundary_conditions" => Any[],
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    @test esm["grids"]["domain"]["family"] == "cartesian"
end

@testitem "GDD family dispatch: latlon family propagates to ESM grid block" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "latlon",
            "spatial" => Dict{String,Any}(
                "lat" => Dict{String,Any}("min" => -90.0, "max" => 90.0,  "grid_spacing" => 5.0),
                "lon" => Dict{String,Any}("min" => 0.0,   "max" => 360.0, "grid_spacing" => 5.0),
            ),
            "boundary_conditions" => Any[
                Dict{String,Any}("type" => "periodic", "dimensions" => ["lon"]),
            ],
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "latlon"
    # lat and lon axes emitted (sorted: lat before lon)
    names = [d["name"] for d in block["dimensions"]]
    @test "lat" in names
    @test "lon" in names
    # periodic flag: lon is periodic, lat is not
    lat_dim = only(filter(d -> d["name"] == "lat", block["dimensions"]))
    lon_dim = only(filter(d -> d["name"] == "lon", block["dimensions"]))
    @test lat_dim["periodic"] == false
    @test lon_dim["periodic"] == true
    @test lat_dim["size"] == 36   # (90 - (-90)) / 5
    @test lon_dim["size"] == 72   # 360 / 5
end

@testitem "GDD family dispatch: vertical family propagates to ESM grid block" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "vertical",
            "spatial" => Dict{String,Any}(
                "k" => Dict{String,Any}("min" => 0.0, "max" => 1.0, "grid_spacing" => 0.1),
            ),
            "boundary_conditions" => Any[],
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "vertical"
    @test block["dimensions"][1]["name"] == "k"
    @test block["dimensions"][1]["size"] == 10
    @test block["dimensions"][1]["periodic"] == false
end

@testitem "GDD family dispatch: arakawa family propagates to ESM grid block" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "arakawa",
            "spatial" => Dict{String,Any}(
                "x" => Dict{String,Any}("min" => 0.0, "max" => 1.0, "grid_spacing" => 0.1),
                "y" => Dict{String,Any}("min" => 0.0, "max" => 1.0, "grid_spacing" => 0.1),
            ),
            "boundary_conditions" => Any[
                Dict{String,Any}("type" => "periodic", "dimensions" => ["x", "y"]),
            ],
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "arakawa"
    names = [d["name"] for d in block["dimensions"]]
    @test "x" in names && "y" in names
end

@testitem "GDD rules dispatch: cartesian stencil → lower_stencil_to_scheme (scheme + use-rule)" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    cartesian_stencil = Dict{String,Any}(
        "centered_2nd_uniform" => Dict{String,Any}(
            "applies_to" => Dict{String,Any}(
                "op" => "grad", "args" => ["\$u"], "dim" => "\$x",
            ),
            "grid_family" => "cartesian",
            "combine"     => "+",
            "stencil"     => Any[
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "cartesian", "axis" => "\$x", "offset" => -1), "coeff" => -0.5),
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "cartesian", "axis" => "\$x", "offset" =>  1), "coeff" =>  0.5),
            ],
        ),
    )
    _inject_rules!(esm, cartesian_stencil, "/dev/null")

    # Cartesian stencil: goes through lower_stencil_to_scheme →
    # a `discretizations` scheme block is emitted + a `use:` rule appended.
    @test haskey(esm["discretizations"], "centered_2nd_uniform")
    scheme = esm["discretizations"]["centered_2nd_uniform"]
    @test scheme["grid_family"] == "cartesian"
    @test length(esm["rules"]) == 1
    @test esm["rules"][1]["use"] == "centered_2nd_uniform"
end

@testitem "GDD rules dispatch: latlon stencil → lower_stencil_to_replacement (pattern+replacement rule)" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    latlon_stencil = Dict{String,Any}(
        "latlon_grad" => Dict{String,Any}(
            "applies_to" => Dict{String,Any}(
                "op" => "grad", "args" => ["\$u"], "dim" => "\$k",
            ),
            "grid_family" => "latlon",
            "combine"     => "+",
            "stencil"     => Any[
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "latlon", "axis" => "lon", "offset" => -1), "coeff" => -0.5),
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "latlon", "axis" => "lon", "offset" =>  1), "coeff" =>  0.5),
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "latlon", "axis" => "lat", "offset" => -1), "coeff" => -0.5),
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "latlon", "axis" => "lat", "offset" =>  1), "coeff" =>  0.5),
            ],
        ),
    )
    _inject_rules!(esm, latlon_stencil, "/dev/null")

    # Latlon stencil: goes through lower_stencil_to_replacement →
    # no scheme block emitted; a plain pattern+replacement rule appended.
    @test !haskey(esm["discretizations"], "latlon_grad")
    @test length(esm["rules"]) == 1
    rule = esm["rules"][1]
    @test rule["name"] == "latlon_grad"
    @test haskey(rule, "replacement")
    @test !haskey(rule, "use")   # must be replacement-form, not scheme-use-form
end

@testitem "GDD rules dispatch: vertical stencil → lower_stencil_to_replacement (pattern+replacement rule)" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    vertical_stencil = Dict{String,Any}(
        "vertical_grad" => Dict{String,Any}(
            "applies_to" => Dict{String,Any}(
                "op" => "grad", "args" => ["\$u"], "dim" => "\$k",
            ),
            "grid_family" => "vertical",
            "combine"     => "+",
            "stencil"     => Any[
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "vertical", "axis" => "\$k", "stagger" => "face_bottom", "offset" => 0), "coeff" => -1.0),
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "vertical", "axis" => "\$k", "stagger" => "face_top",    "offset" => 0), "coeff" =>  1.0),
            ],
        ),
    )
    _inject_rules!(esm, vertical_stencil, "/dev/null")

    # Vertical stencil: goes through lower_stencil_to_replacement →
    # no scheme block emitted; a plain pattern+replacement rule appended.
    @test !haskey(esm["discretizations"], "vertical_grad")
    @test length(esm["rules"]) == 1
    rule = esm["rules"][1]
    @test rule["name"] == "vertical_grad"
    @test haskey(rule, "replacement")
    @test !haskey(rule, "use")
end

@testitem "GDD grids dispatch: mpas family emits n_cells dimension" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "mpas",
            "n_cells" => 642,
            "loader"  => Dict{String,Any}("path" => "meshes/x1.642.grid.nc", "reader" => "mpas_mesh"),
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "mpas"
    cell_dim = only(filter(d -> d["name"] == "n_cells", block["dimensions"]))
    @test cell_dim["size"] == 642
    @test cell_dim["periodic"] == false
    @test cell_dim["spacing"] == "unstructured"
end

@testitem "GDD grids dispatch: mpas family with n_edges emits both dimensions" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "mpas",
            "n_cells" => 642,
            "n_edges" => 1920,
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "mpas"
    names = [d["name"] for d in block["dimensions"]]
    @test "n_cells" in names
    @test "n_edges" in names
    cell_dim = only(filter(d -> d["name"] == "n_cells", block["dimensions"]))
    @test cell_dim["size"] == 642
end

@testitem "GDD grids dispatch: mpas missing n_cells raises" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}("family" => "mpas"),
    )
    @test_throws ArgumentError _inject_grids!(esm, gdd_grids, "/dev/null")
end

@testitem "GDD grids dispatch: duo family emits n_cells dimension" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "family"  => "duo",
            "n_cells" => 320,
            "loader"  => Dict{String,Any}("path" => "builtin://icosahedral/2", "reader" => "auto"),
        ),
    )
    _inject_grids!(esm, gdd_grids, "/dev/null")

    block = esm["grids"]["domain"]
    @test block["family"] == "duo"
    cell_dim = only(filter(d -> d["name"] == "n_cells", block["dimensions"]))
    @test cell_dim["size"] == 320
    @test cell_dim["periodic"] == false
    @test cell_dim["spacing"] == "unstructured"
end

@testitem "GDD grids dispatch: duo family mesh ladder (levels 2/3/4)" begin
    using EarthSciDiscretizations: _inject_grids!

    for (level, expected_n) in ((2, 320), (3, 1280), (4, 5120))
        esm = Dict{String,Any}("grids" => Dict{String,Any}())
        gdd_grids = Dict{String,Any}(
            "domain" => Dict{String,Any}(
                "family"  => "duo",
                "n_cells" => expected_n,
                "loader"  => Dict{String,Any}("path" => "builtin://icosahedral/$level", "reader" => "auto"),
            ),
        )
        _inject_grids!(esm, gdd_grids, "/dev/null")
        block = esm["grids"]["domain"]
        @test block["family"] == "duo"
        cell_dim = only(filter(d -> d["name"] == "n_cells", block["dimensions"]))
        @test cell_dim["size"] == expected_n
    end
end

@testitem "GDD grids dispatch: duo missing n_cells raises" begin
    using EarthSciDiscretizations: _inject_grids!

    esm = Dict{String,Any}("grids" => Dict{String,Any}())
    gdd_grids = Dict{String,Any}(
        "domain" => Dict{String,Any}("family" => "duo"),
    )
    err = @test_throws ArgumentError _inject_grids!(esm, gdd_grids, "/dev/null")
    @test occursin("DUO GDD", sprint(showerror, err.value))
    @test occursin("n_cells", sprint(showerror, err.value))
end

@testitem "GDD rules dispatch: reduction selector throws ESS/esm-bpr gate error" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    unstructured_stencil = Dict{String,Any}(
        "nn_diffusion_mpas" => Dict{String,Any}(
            "applies_to"  => Dict{String,Any}("op" => "laplacian", "args" => ["\$u"], "dim" => "cell"),
            "grid_family" => "unstructured",
            "combine"     => "+",
            "stencil"     => Any[
                Dict{String,Any}("selector" => Dict{String,Any}("kind" => "reduction", "table" => "cells_on_cell", "count_expr" => 6, "k_bound" => "k", "combine" => "+"), "coeff" => 1.0),
            ],
        ),
    )
    err = @test_throws ArgumentError _inject_rules!(esm, unstructured_stencil, "/dev/null")
    @test occursin("ESS/esm-bpr", sprint(showerror, err.value))
end
