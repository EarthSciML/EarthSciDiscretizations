using Test
using TestItems

@testitem "lower_stencil_to_replacement: passes replacement-form rules through" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
        "replacement" => Dict("op" => "+", "args" => Any[1, 2]),
    )
    out = lower_stencil_to_replacement(rule)
    @test out["replacement"] == rule["replacement"]
    # idempotent: re-running yields the same replacement
    @test lower_stencil_to_replacement(out)["replacement"] == rule["replacement"]
end


@testitem "lower_stencil_to_replacement: ESS parse_expression accepts lowered upwind_1st AST" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # The lowered `replacement` must be a well-formed ExpressionNode AST that
    # ESS's `parse_expression` accepts. This guards the lowerer's output
    # contract — if the structure drifts, ESS rule_engine will reject the
    # rule at `parse_rule` time, so catching it at parse_expression here is
    # the cheapest signal that the AST shape is right.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st"])
    out = lower_stencil_to_replacement(rule)

    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing
end

# esd-7h2: removed "lower_stencil_to_replacement: lowers cartesian upwind_1st_nonuniform (esd-02z)"
# upwind_1st_nonuniform.json was migrated to arrayop replacement form — it no longer has a
# stencil array, so the stencil-lowering path is not exercised by this rule.
# The rule is tested via Layer-A (canonical byte contract) and Layer-B (MMS SKIP, no mms_kind).

# esd-7h2: removed "lower_stencil_to_replacement: ESS parse_expression accepts lowered upwind_1st_nonuniform AST (esd-02z)"
# Same migration as above — rule is now in replacement form, stencil path retired.

@testitem "lower_stencil_to_replacement: errors on missing stencil and replacement" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
    )
    @test_throws ArgumentError lower_stencil_to_replacement(rule)
end

@testitem "lower_stencil_to_replacement: errors on missing applies_to" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "cartesian", "axis" => "\$x", "offset" => 0),
            "coeff" => 1,
        )],
    )
    @test_throws ArgumentError lower_stencil_to_replacement(rule)
end

@testitem "lower_stencil_to_replacement: errors on unsupported selector kind" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "indirect", "axis" => "n", "offset" => 0),
            "coeff" => 1,
        )],
    )
    err = try
        lower_stencil_to_replacement(rule)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("'indirect'", err.msg)
    @test occursin("not yet supported", err.msg)
end

# esd-7h2: removed "lower_stencil_to_replacement: errors on axis mismatch"
# The axis-mismatch check was inside _lower_cartesian_entry (retired with the
# cartesian branch). Vertical and latlon kinds validate their own axis fields
# independently; no replacement test needed for the retired cartesian path.

@testitem "lower_stencil_to_replacement: single entry produces bare term (vertical)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    # Uses vertical stencil (cell_center stagger, offset 1) to verify the
    # general "single entry → bare term" behaviour (no combine wrapper).
    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "vertical", "axis" => "\$k",
                               "stagger" => "cell_center", "offset" => 1),
            "coeff" => 7,
        )],
    )
    out = lower_stencil_to_replacement(rule)
    repl = out["replacement"]
    # No combine wrapper for a single entry
    @test repl["op"] == "*"
    @test repl["args"][1] == 7
    idx = repl["args"][2]
    @test idx["op"] == "index"
    @test idx["args"][2]["op"] == "+"
    @test idx["args"][2]["args"][2] == 1
end

@testitem "lower_stencil_to_replacement: lowers divergence_arakawa_c" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    # `divergence_arakawa_c` carries a canonical bare replacement (esd-eg5): no stencil,
    # no arrayop wrapper. The replacement is a "+" of 4 coeff*index terms using i/j axes.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "divergence_arakawa_c.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["divergence_arakawa_c"])

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 4

    # All 4 terms are coeff * index(field, i_arg, j_arg)
    for term in repl["args"]
        @test term["op"] == "*"
        idx = term["args"][2]
        @test idx["op"] == "index"
        @test length(idx["args"]) == 3
    end

    # Entry 1: index($Fx, i, j), coeff = -1/dx
    e1 = repl["args"][1]
    idx1 = e1["args"][2]
    @test String(idx1["args"][1]) == "\$Fx"
    @test String(idx1["args"][2]) == "i"
    @test String(idx1["args"][3]) == "j"

    # Entry 2: index($Fx, i+1, j), coeff = +1/dx
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test String(idx2["args"][1]) == "\$Fx"
    a2_i = idx2["args"][2]
    @test a2_i["op"] == "+"
    @test String(a2_i["args"][1]) == "i"
    @test Int(a2_i["args"][2]) == 1
    @test String(idx2["args"][3]) == "j"

    # Entry 3: index($Fy, i, j), coeff = -1/dy
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    @test String(idx3["args"][1]) == "\$Fy"
    @test String(idx3["args"][2]) == "i"
    @test String(idx3["args"][3]) == "j"

    # Entry 4: index($Fy, i, j+1), coeff = +1/dy
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    @test String(idx4["args"][1]) == "\$Fy"
    @test String(idx4["args"][2]) == "i"
    a4_j = idx4["args"][3]
    @test a4_j["op"] == "+"
    @test String(a4_j["args"][1]) == "j"
    @test Int(a4_j["args"][2]) == 1

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: ESS parse_expression accepts lowered divergence_arakawa_c" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "divergence_arakawa_c.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["divergence_arakawa_c"])
    out = lower_stencil_to_replacement(rule)

    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing
end

@testitem "lower_stencil_to_replacement: errors on mixed selector kinds" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[
            Dict(
                "selector" => Dict(
                    "kind" => "arakawa", "stagger" => "face_x",
                    "axis" => "\$x", "offset" => 0,
                ),
                "coeff" => 1,
            ),
            Dict(
                "selector" => Dict("kind" => "cartesian", "axis" => "\$x", "offset" => 0),
                "coeff" => 1,
            ),
        ],
    )
    err = try
        lower_stencil_to_replacement(rule)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("mixes selector kinds", err.msg)
end

# esd-7h2: removed "lower_stencil_to_replacement: lowers lax_friedrichs_flux (cartesian)"
# lax_friedrichs_flux.json retains its cartesian stencil (EINSUM-8 scope) but
# _lower_cartesian_entry is retired; the rule is FV and does not flow through
# lower_stencil_to_replacement in any active code path. Tested via Layer-A skip
# (applicable:false in convergence fixture) and the hand-pinned flux test in
# test_lax_friedrichs_flux_rule.jl.

# esd-3d7: removed "lower_stencil_to_replacement: lowers centered_2nd_deriv_uniform (cartesian)"
# centered_2nd_deriv_uniform.json was migrated to arrayop replacement form — it no longer has a
# stencil array, so the stencil-lowering path is not exercised by this rule.
# The rule is tested via Layer-A (canonical byte contract) and Layer-B (MMS convergence).

@testitem "lower_stencil_to_replacement: lowers laplacian_2nd_uniform_cartesian (arakawa)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # laplacian_2nd_uniform_cartesian carries a canonical bare replacement (esd-eg5):
    # no stencil, no arrayop wrapper. Six coeff*index terms using i/j axes.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "laplacian_2nd_uniform_cartesian.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["laplacian_2nd_uniform_cartesian"])
    @test rule["applies_to"]["op"] == "laplacian"
    @test String(rule["applies_to"]["args"][1]) == "\$u"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 6

    # All terms are coeff * index($u, i_arg, j_arg)
    for term in repl["args"]
        @test term["op"] == "*"
        idx = term["args"][2]
        @test idx["op"] == "index"
        @test length(idx["args"]) == 3
        @test String(idx["args"][1]) == "\$u"
    end

    # Entries 1–3: shift i axis (x direction)
    # Entry 1: index($u, i+(-1), j)
    e1 = repl["args"][1]
    idx1 = e1["args"][2]
    @test idx1["args"][2]["op"] == "+"
    @test String(idx1["args"][2]["args"][1]) == "i"
    @test Int(idx1["args"][2]["args"][2]) == -1
    @test String(idx1["args"][3]) == "j"

    # Entry 2: index($u, i, j) — bare i
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test String(idx2["args"][2]) == "i"
    @test String(idx2["args"][3]) == "j"

    # Entry 3: index($u, i+1, j)
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    @test idx3["args"][2]["op"] == "+"
    @test Int(idx3["args"][2]["args"][2]) == 1
    @test String(idx3["args"][3]) == "j"

    # Entries 4–6: shift j axis (y direction)
    # Entry 4: index($u, i, j+(-1))
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    @test String(idx4["args"][2]) == "i"
    @test idx4["args"][3]["op"] == "+"
    @test String(idx4["args"][3]["args"][1]) == "j"
    @test Int(idx4["args"][3]["args"][2]) == -1

    # Entry 5: index($u, i, j) — bare j
    e5 = repl["args"][5]
    idx5 = e5["args"][2]
    @test String(idx5["args"][2]) == "i"
    @test String(idx5["args"][3]) == "j"

    # Entry 6: index($u, i, j+1)
    e6 = repl["args"][6]
    idx6 = e6["args"][2]
    @test String(idx6["args"][2]) == "i"
    @test idx6["args"][3]["op"] == "+"
    @test Int(idx6["args"][3]["args"][2]) == 1

    # ESS parse_expression accepts the replacement AST directly
    ess_expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test ess_expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: lowers staggered_1st_uniform_cc_to_face" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # staggered_1st_uniform_cc_to_face carries an explicit arrayop replacement AST since
    # EINSUM-4 (esd-770). 1D x-axis, offsets {-1, 0}. Mirrors MOL StaggeredStencilInfo
    # EdgeAlignedVar: d/dx(u)|_face[j] = (u_cc[j] - u_cc[j-1]) / h.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "staggered_1st_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["staggered_1st_uniform_cc_to_face"])

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "arrayop"
    @test repl["output_idx"] == ["\$x"]
    expr = repl["expr"]
    @test expr["op"] == "+"
    @test length(expr["args"]) == 2

    # Entry 1: cell_center, axis $x, offset -1 -> index($u, $x + (-1)), coeff -1/h
    e1 = expr["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    arg1 = idx1["args"][2]
    @test arg1["op"] == "+"
    @test String(arg1["args"][1]) == "\$x"
    @test Int(arg1["args"][2]) == -1

    # Entry 2: cell_center, axis $x, offset 0 -> index($u, $x) (no `+ 0` wrapper)
    e2 = expr["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$x"

    # ESS parse_expression accepts the inner expr AST
    ess_expr = EarthSciSerialization.parse_expression(expr)
    @test ess_expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: lowers staggered_1st_uniform_face_to_cc" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # staggered_1st_uniform_face_to_cc carries an explicit arrayop replacement AST since
    # EINSUM-4 (esd-770). 1D x-axis, offsets {0, +1}. Mirrors MOL StaggeredStencilInfo
    # CenterAlignedVar: d/dx(u)|_cc[i] = (u_face[i+1] - u_face[i]) / h.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "staggered_1st_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["staggered_1st_uniform_face_to_cc"])

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "arrayop"
    @test repl["output_idx"] == ["\$x"]
    expr = repl["expr"]
    @test expr["op"] == "+"
    @test length(expr["args"]) == 2

    # Entry 1: face_x, axis $x, offset 0 -> index($u, $x) (no `+ 0` wrapper), coeff -1/h
    e1 = expr["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    @test String(idx1["args"][2]) == "\$x"

    # Entry 2: face_x, axis $x, offset 1 -> index($u, $x + 1), coeff +1/h
    e2 = expr["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    arg2 = idx2["args"][2]
    @test arg2["op"] == "+"
    @test String(arg2["args"][1]) == "\$x"
    @test Int(arg2["args"][2]) == 1

    # ESS parse_expression accepts the inner expr AST
    ess_expr = EarthSciSerialization.parse_expression(expr)
    @test ess_expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

# =============================================================================
# Latlon family tests
# =============================================================================

@testitem "lower_stencil_to_replacement: lowers centered_2nd_uniform_latlon (latlon)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "centered_2nd_uniform_latlon.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["centered_2nd_uniform_latlon"])
    @test !haskey(rule, "replacement")
    @test rule["grid_family"] == "latlon"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    # 4 stencil entries (2 for lon, 2 for lat) → combine(4 terms)
    @test repl["op"] == "+"
    @test length(repl["args"]) == 4

    # All terms: coeff * index($u, lat_arg, lon_arg) — axes sorted: lat < lon
    for (i, term) in enumerate(repl["args"])
        @test term["op"] == "*"
        idx = term["args"][2]
        @test idx["op"] == "index"
        @test length(idx["args"]) == 3  # operand + lat_arg + lon_arg
        @test String(idx["args"][1]) == "\$u"
    end

    # Entry 1: axis=lon, offset=-1 → index($u, "lat", "lon" + (-1))
    e1 = repl["args"][1]
    idx1 = e1["args"][2]
    @test String(idx1["args"][2]) == "lat"   # lat bare (not the shifted axis)
    a1_lon = idx1["args"][3]
    @test a1_lon["op"] == "+"
    @test String(a1_lon["args"][1]) == "lon"
    @test Int(a1_lon["args"][2]) == -1

    # Entry 2: axis=lon, offset=+1 → index($u, "lat", "lon" + 1)
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test String(idx2["args"][2]) == "lat"
    a2_lon = idx2["args"][3]
    @test a2_lon["op"] == "+"
    @test String(a2_lon["args"][1]) == "lon"
    @test Int(a2_lon["args"][2]) == 1

    # Entry 3: axis=lat, offset=-1 → index($u, "lat" + (-1), "lon")
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    a3_lat = idx3["args"][2]
    @test a3_lat["op"] == "+"
    @test String(a3_lat["args"][1]) == "lat"
    @test Int(a3_lat["args"][2]) == -1
    @test String(idx3["args"][3]) == "lon"   # lon bare

    # Entry 4: axis=lat, offset=+1 → index($u, "lat" + 1, "lon")
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    a4_lat = idx4["args"][2]
    @test a4_lat["op"] == "+"
    @test String(a4_lat["args"][1]) == "lat"
    @test Int(a4_lat["args"][2]) == 1
    @test String(idx4["args"][3]) == "lon"

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: latlon errors on \$-prefixed axis" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "latlon", "axis" => "\$lon", "offset" => 0),
            "coeff" => 1,
        )],
    )
    err = try
        lower_stencil_to_replacement(rule)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("literal geographic axis name", err.msg)
end

# =============================================================================
# Cubed-sphere family tests
# =============================================================================

@testitem "covariant_laplacian_cubed_sphere: cross_metric schema + per-axis stencil lowering (esd-p81)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "covariant_laplacian_cubed_sphere.json",
    )
    raw = JSON.parsefile(path)
    disc = raw["discretizations"]

    # Top-level rule is now a cross_metric composite — not a stencil rule.
    rule = Dict{String, Any}(disc["covariant_laplacian_cubed_sphere"])
    @test rule["kind"] == "cross_metric"
    @test rule["grid_family"] == "cubed_sphere"
    @test !haskey(rule, "stencil")
    @test !haskey(rule, "replacement")

    # 8 cross_metric terms: 2 diagonal second-deriv, 2 cross-deriv, 4 metric-gradient corrections.
    terms = rule["terms"]
    @test length(terms) == 8

    # Term metric_component names match the full covariant Laplacian expansion.
    mc = [t["metric_component"] for t in terms]
    @test "ginv_xi_xi"  in mc
    @test "ginv_eta_eta" in mc
    @test "ginv_xi_eta"  in mc   # two entries (ξη and ηξ for symmetric metric)
    @test "dJgxx_dxi"   in mc
    @test "dJgxe_deta"  in mc
    @test "dJgyy_deta"  in mc
    @test "dJgxe_dxi"   in mc

    # All axis_stencil references resolve to entries in the same discretizations block.
    expected_per_axis = Set([
        "d2_dxi2_cubed_sphere",
        "d2_deta2_cubed_sphere",
        "d2_dxieta_cubed_sphere",
        "d1_dxi_over_J_cubed_sphere",
        "d1_deta_over_J_cubed_sphere",
    ])
    referenced = Set(t["axis_stencil"] for t in terms)
    @test referenced ⊆ expected_per_axis
    for name in referenced
        @test haskey(disc, name)
    end

    # Per-axis stencils are still in stencil form and can be lowered to replacement.
    for name in expected_per_axis
        pa = Dict{String, Any}(disc[name])
        @test haskey(pa, "stencil")
        @test !haskey(pa, "replacement")
        out = lower_stencil_to_replacement(pa)
        @test haskey(out, "replacement")
        expr = EarthSciSerialization.parse_expression(out["replacement"])
        @test expr !== nothing
    end

    # lower_stencil_to_replacement throws on the cross_metric top-level rule
    # (it has neither stencil nor replacement).
    @test_throws ArgumentError lower_stencil_to_replacement(rule)
end

@testitem "lower_stencil_to_replacement: lowers transport_2d (cubed_sphere, finite_volume)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "transport_2d.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["transport_2d"])
    @test !haskey(rule, "replacement")
    @test rule["grid_family"] == "cubed_sphere"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    # 5 stencil entries → combine(5 terms); axes: eta < xi
    @test repl["op"] == "+"
    @test length(repl["args"]) == 5

    for term in repl["args"]
        @test term["op"] == "*"
        idx = term["args"][2]
        @test idx["op"] == "index"
        @test length(idx["args"]) == 3
        @test String(idx["args"][1]) == "\$q"
    end

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

# =============================================================================
# Vertical family tests
# =============================================================================

@testitem "lower_stencil_to_replacement: lowers centered_2nd_uniform_vertical (vertical)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "centered_2nd_uniform_vertical.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["centered_2nd_uniform_vertical"])
    @test !haskey(rule, "replacement")
    @test rule["grid_family"] == "vertical"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    # 2 stencil entries → combine(2 terms)
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # Both entries: coeff * index($u, $k) — offset=0, no `+ 0` wrapper
    # Stagger (face_bottom / face_top) is NOT in the AST
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    @test String(idx1["args"][2]) == "\$k"   # offset=0 → bare axis var

    e2 = repl["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$k"   # offset=0 → bare axis var

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: vertical errors on bad stagger" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
        "stencil" => Any[Dict(
            "selector" => Dict(
                "kind" => "vertical", "axis" => "\$k",
                "stagger" => "unknown_stagger", "offset" => 0,
            ),
            "coeff" => 1,
        )],
    )
    err = try
        lower_stencil_to_replacement(rule)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("'unknown_stagger'", err.msg)
end

@testitem "lower_stencil_to_replacement: vertical errors on axis mismatch" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
        "stencil" => Any[Dict(
            "selector" => Dict(
                "kind" => "vertical", "axis" => "\$z",
                "stagger" => "face_bottom", "offset" => 0,
            ),
            "coeff" => 1,
        )],
    )
    err = try
        lower_stencil_to_replacement(rule)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("disagrees", err.msg)
end

# esd-3d7: removed "lower_stencil_to_scheme: lowers cartesian upwind_1st to ESS scheme parts" and
# "lower_stencil_to_scheme: ESS parse_schemes + parse_rule accept the lowered parts".
# upwind_1st.json was migrated to arrayop replacement form — it no longer has a stencil array, so
# lower_stencil_to_scheme throws ArgumentError on it. The rule is exercised via Layer-A and Layer-B.

@testitem "lower_stencil_to_scheme: errors on non-cartesian grid_family and selector kind" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme

    base = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
        "grid_family" => "cartesian",
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "cartesian", "axis" => "\$x", "offset" => 1),
            "coeff" => 1,
        )],
    )

    vertical = copy(base); vertical["grid_family"] = "vertical"
    err = try lower_stencil_to_scheme("r", vertical); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("not lowerable", err.msg)

    bad_kind = copy(base)
    bad_kind["stencil"] = Any[Dict(
        "selector" => Dict("kind" => "arakawa", "axis" => "\$x", "offset" => 1,
                           "stagger" => "cell_center"),
        "coeff" => 1,
    )]
    err = try lower_stencil_to_scheme("r", bad_kind); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("'arakawa'", err.msg)

    no_stencil = Dict{String, Any}(
        "applies_to" => base["applies_to"], "grid_family" => "cartesian",
        "replacement" => Dict("op" => "+", "args" => Any[1, 1]),
    )
    err = try lower_stencil_to_scheme("r", no_stencil); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("no 'stencil'", err.msg)

    axis_mismatch = copy(base)
    axis_mismatch["stencil"] = Any[Dict(
        "selector" => Dict("kind" => "cartesian", "axis" => "\$y", "offset" => 1),
        "coeff" => 1,
    )]
    err = try lower_stencil_to_scheme("r", axis_mismatch); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("disagrees", err.msg)
end

@testitem "lower_stencil_to_scheme: lowers unstructured nn_diffusion_mpas to ESS scheme parts" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "nn_diffusion_mpas.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["nn_diffusion_mpas"])

    scheme, use_rule = lower_stencil_to_scheme("nn_diffusion_mpas", rule)

    @test scheme["grid_family"] == "unstructured"
    @test scheme["combine"] == "+"
    @test scheme["applies_to"] == rule["applies_to"]
    @test length(scheme["stencil"]) == 2

    sel1 = scheme["stencil"][1]["selector"]
    @test sel1["kind"] == "reduction"
    @test sel1["table"] == "cells_on_cell"
    @test sel1["k_bound"] == "k"
    @test sel1["combine"] == "+"

    sel2 = scheme["stencil"][2]["selector"]
    @test sel2["kind"] == "indirect"
    @test sel2["table"] == "cells_on_cell"
    @test sel2["index_expr"] == "\$target"

    @test use_rule == Dict{String, Any}(
        "name"    => "nn_diffusion_mpas",
        "pattern" => rule["applies_to"],
        "use"     => "nn_diffusion_mpas",
    )

    # ESS parse_scheme must accept the emitted scheme without error.
    schemes = EarthSciSerialization.parse_schemes(Dict("nn_diffusion_mpas" => scheme))
    @test haskey(schemes, "nn_diffusion_mpas")
    sch = schemes["nn_diffusion_mpas"]
    @test sch.grid_family == "unstructured"
    @test length(sch.stencil) == 2
end

@testitem "lower_stencil_to_scheme: unstructured discretize() produces arrayop equations" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "nn_diffusion_mpas.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["nn_diffusion_mpas"])
    scheme, use_rule = lower_stencil_to_scheme("nn_diffusion_mpas", rule)

    esm = Dict{String, Any}(
        "esm"   => "0.5.0",
        "grids" => Dict{String, Any}(
            "domain" => Dict{String, Any}(
                "family"     => "unstructured",
                "dimensions" => Any[Dict{String, Any}("name" => "n_cells", "size" => 5)],
            ),
        ),
        "models" => Dict{String, Any}(
            "diffusion" => Dict{String, Any}(
                "grid" => "domain",
                "variables" => Dict{String, Any}(
                    "u" => Dict{String, Any}(
                        "type"     => "state",
                        "default"  => 0.0,
                        "units"    => "1",
                        "shape"    => Any["n_cells"],
                        "location" => "cell_center",
                    ),
                ),
                "equations" => Any[
                    Dict{String, Any}(
                        "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u"], "wrt" => "t"),
                        "rhs" => Dict{String, Any}("op" => "laplacian", "args" => Any["u"], "dim" => "cell"),
                    ),
                ],
            ),
        ),
        "discretizations" => Dict{String, Any}("nn_diffusion_mpas" => scheme),
        "rules" => Any[use_rule],
    )

    out = EarthSciSerialization.discretize(esm; strict_unrewritten = true)
    eqs = out["models"]["diffusion"]["equations"]
    @test !isempty(eqs)
    rhs = eqs[1]["rhs"]
    @test rhs isa AbstractDict
    # The laplacian was rewritten (rule applied); RHS is the expanded combine of
    # reduction + indirect terms — not the bare "laplacian" op any more.
    @test rhs["op"] != "laplacian"
    # At least one arrayop node appears in the expanded RHS.
    @test occursin("arrayop", JSON.json(rhs))
end

@testitem "lower_stencil_to_scheme: multi-output stencil object lowers per output" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "ppm_reconstruction.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["ppm_reconstruction"])

    # Without output= the multi-output object is ambiguous.
    err = try lower_stencil_to_scheme("ppm_reconstruction", rule); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("multi-output", err.msg)
    @test occursin("q_left_edge", err.msg)

    # Unknown output names are rejected with the available set.
    err = try lower_stencil_to_scheme("ppm_reconstruction", rule; output = "nope"); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("q_right_edge", err.msg)

    # Each named output lowers to an ordinary single-axis cartesian scheme.
    scheme, use_rule = lower_stencil_to_scheme("ppm_reconstruction", rule; output = "q_left_edge")
    @test scheme["grid_family"] == "cartesian"
    @test length(scheme["stencil"]) == 4
    offsets = [e["selector"]["offset"] for e in scheme["stencil"]]
    @test offsets == [-2, -1, 0, 1]
    @test use_rule["pattern"] == rule["applies_to"]
    @test use_rule["use"] == "ppm_reconstruction"

    scheme_r, _ = lower_stencil_to_scheme("ppm_reconstruction", rule; output = "q_right_edge")
    @test [e["selector"]["offset"] for e in scheme_r["stencil"]] == [-1, 0, 1, 2]

    # output= on a flat single-output stencil is an authoring error.
    flat = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
        "grid_family" => "cartesian",
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "cartesian", "axis" => "\$x", "offset" => 0),
            "coeff" => 1,
        )],
    )
    err = try lower_stencil_to_scheme("r", flat; output = "q_left_edge"); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("flat", err.msg)
end

