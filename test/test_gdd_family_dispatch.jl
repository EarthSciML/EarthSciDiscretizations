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

# esd-7h2: removed "GDD rules dispatch: cartesian stencil → lower_stencil_to_scheme (scheme + use-rule)"
# The cartesian dispatch branch was retired: upwind_1st_nonuniform.json and
# centered_2nd_nonuniform_cartesian.json were migrated to arrayop replacement form,
# and the `if kind == "cartesian"` branch in _inject_rules! was removed.
# All remaining FD catalog rules with cartesian selectors are already in replacement
# or arrayop form and take the `elseif haskey(spec, "replacement")` path instead.

# esd-t4h: removed "GDD rules dispatch: latlon stencil → lower_stencil_to_replacement (pattern+replacement rule)".
# The latlon stencil dispatch branch was retired in esd-t4h: centered_2nd_uniform_latlon.json now
# carries an authored replacement AST (no stencil), and the stencil dispatch for kind="latlon" was
# removed from _inject_rules!. The latlon authored-replacement path is tested below.

@testitem "GDD rules dispatch: latlon authored replacement → axis substitution (lat→i, lon→j)" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    latlon_replacement = Dict{String,Any}(
        "latlon_grad" => Dict{String,Any}(
            "applies_to"  => Dict{String,Any}("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
            "grid_family" => "latlon",
            "replacement" => Dict{String,Any}(
                "op"   => "+",
                "args" => Any[
                    Dict{String,Any}("op" => "*", "args" => Any[-0.5,
                        Dict{String,Any}("op" => "index", "args" => Any["\$u", "lat",
                            Dict{String,Any}("op" => "+", "args" => Any["lon", -1])])]),
                    Dict{String,Any}("op" => "*", "args" => Any[ 0.5,
                        Dict{String,Any}("op" => "index", "args" => Any["\$u", "lat",
                            Dict{String,Any}("op" => "+", "args" => Any["lon",  1])])]),
                ],
            ),
        ),
    )
    _inject_rules!(esm, latlon_replacement, "/dev/null")

    # Authored replacement: no scheme block; a pattern+replacement rule appended.
    @test !haskey(esm["discretizations"], "latlon_grad")
    @test length(esm["rules"]) == 1
    rule = esm["rules"][1]
    @test rule["name"] == "latlon_grad"
    @test haskey(rule, "replacement")
    @test !haskey(rule, "use")

    # Geographic axis names ("lat", "lon") must be translated to canonical
    # loop-variable names ("i", "j") so build_evaluator resolves index expressions.
    repl = rule["replacement"]
    idx1 = repl["args"][1]["args"][2]  # first term's index node
    @test String(idx1["args"][2]) == "i"   # lat → i
    lon_arg = idx1["args"][3]
    @test lon_arg["op"] == "+"
    @test String(lon_arg["args"][1]) == "j"   # lon → j
end

# esd-t4h: removed "GDD rules dispatch: vertical stencil → lower_stencil_to_replacement (pattern+replacement rule)".
# The vertical stencil dispatch branch was retired in esd-t4h: centered_2nd_uniform_vertical.json now
# carries an authored replacement AST (no stencil), and the stencil dispatch for kind="vertical" was
# removed from _inject_rules!. The vertical authored-replacement path is tested below.

@testitem "GDD rules dispatch: vertical authored replacement → replacement rule, no axis subst" begin
    using EarthSciDiscretizations: _inject_rules!

    esm = Dict{String,Any}("rules" => Any[], "discretizations" => Dict{String,Any}())
    vertical_replacement = Dict{String,Any}(
        "vertical_grad" => Dict{String,Any}(
            "applies_to"  => Dict{String,Any}("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
            "grid_family" => "vertical",
            "replacement" => Dict{String,Any}(
                "op"   => "/",
                "args" => Any[
                    Dict{String,Any}("op" => "+", "args" => Any[
                        Dict{String,Any}("op" => "index", "args" => Any["\$u", Dict{String,Any}("op" => "+", "args" => Any["\$k", 1])]),
                        Dict{String,Any}("op" => "-", "args" => Any[
                            Dict{String,Any}("op" => "index", "args" => Any["\$u", Dict{String,Any}("op" => "-", "args" => Any["\$k", 1])]),
                        ]),
                    ]),
                    Dict{String,Any}("op" => "*", "args" => Any[2, "h"]),
                ],
            ),
        ),
    )
    _inject_rules!(esm, vertical_replacement, "/dev/null")

    # Authored replacement: no scheme block; a pattern+replacement rule appended.
    @test !haskey(esm["discretizations"], "vertical_grad")
    @test length(esm["rules"]) == 1
    rule = esm["rules"][1]
    @test rule["name"] == "vertical_grad"
    @test haskey(rule, "replacement")
    @test !haskey(rule, "use")

    # Vertical replacement uses pattern variables ($k); no axis substitution applied.
    repl = rule["replacement"]
    @test repl["op"] == "/"
    num = repl["args"][1]
    @test num["op"] == "+"
    idx_plus = num["args"][1]
    @test idx_plus["op"] == "index"
    kplus = idx_plus["args"][2]
    @test kplus["op"] == "+"
    @test String(kplus["args"][1]) == "\$k"
    @test Int(kplus["args"][2]) == 1
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

# esd-agh: removed "GDD rules dispatch: reduction selector lowers to ESS scheme (ess-t0z)"
# The unstructured stencil dispatch branch (kind="reduction"/"indirect") was retired in esd-agh:
# nn_diffusion_mpas.json and nn_diffusion_duo.json were migrated to authored reduce-arrayop
# replacement form, and lower_stencil_to_scheme + its dispatch in _inject_rules! were deleted.
# Rules with a "replacement" field take the replacement path; stencil-only rules are now silently
# skipped. The unstructured reduction path is tested end-to-end via build_ode_problem in
# test_build_ode_problem.jl (MPAS analytic Laplacian check, esd-aij/esd-agh).
