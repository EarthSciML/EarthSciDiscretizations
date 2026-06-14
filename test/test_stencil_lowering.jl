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

@testitem "lower_stencil_to_replacement: lowers cartesian upwind_1st" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st"])
    # Strip the replacement field (added by esd-0ip for the Layer-B runner) so
    # this test exercises the actual stencil-lowering path and verifies the AST shape.
    stencil_only = Dict{String, Any}(k => v for (k, v) in rule if k != "replacement")
    @test haskey(stencil_only, "stencil")
    @test !haskey(stencil_only, "replacement")

    out = lower_stencil_to_replacement(stencil_only)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # Entry 1: offset = -1, coeff = -1/dx -> -1/dx * index($u, $x + (-1))
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    @test length(e1["args"]) == 2
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    arg1 = idx1["args"][2]
    @test arg1["op"] == "+"
    @test String(arg1["args"][1]) == "\$x"
    @test Int(arg1["args"][2]) == -1

    # Entry 2: offset = 0 -> index($u, $x) (no `+ 0` wrapper)
    e2 = repl["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$x"
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

@testitem "lower_stencil_to_replacement: lowers cartesian upwind_1st_nonuniform (esd-02z)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st_nonuniform.json",
    )
    raw  = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st_nonuniform"])
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # Entry 1: offset = -1, coeff = -1/index(dx,$x) -> term * index($u, $x+(-1))
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    coeff1 = e1["args"][1]
    @test coeff1["op"] == "/"
    # numerator is -1
    @test Int(coeff1["args"][1]) == -1
    # denominator is index("dx", "$x")
    denom1 = coeff1["args"][2]
    @test denom1["op"] == "index"
    @test String(denom1["args"][1]) == "dx"
    @test String(denom1["args"][2]) == "\$x"
    # operand index
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    arg1 = idx1["args"][2]
    @test arg1["op"] == "+"
    @test String(arg1["args"][1]) == "\$x"
    @test Int(arg1["args"][2]) == -1

    # Entry 2: offset = 0, coeff = 1/index(dx,$x) -> term * index($u, $x)
    e2 = repl["args"][2]
    @test e2["op"] == "*"
    coeff2 = e2["args"][1]
    @test coeff2["op"] == "/"
    @test Int(coeff2["args"][1]) == 1
    denom2 = coeff2["args"][2]
    @test denom2["op"] == "index"
    @test String(denom2["args"][1]) == "dx"
    @test String(denom2["args"][2]) == "\$x"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$x"
end

@testitem "lower_stencil_to_replacement: ESS parse_expression accepts lowered upwind_1st_nonuniform AST (esd-02z)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st_nonuniform.json",
    )
    raw  = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st_nonuniform"])
    out  = lower_stencil_to_replacement(rule)

    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing
end

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

@testitem "lower_stencil_to_replacement: errors on axis mismatch" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "cartesian", "axis" => "\$y", "offset" => 1),
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

@testitem "lower_stencil_to_replacement: single entry produces bare term" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$x"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "cartesian", "axis" => "\$x", "offset" => 2),
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
    @test idx["args"][2]["args"][2] == 2
end

@testitem "lower_stencil_to_replacement: lowers divergence_arakawa_c" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    # `divergence_arakawa_c` uses two args: `\$Fx` (face_x component) and
    # `\$Fy` (face_y component), per SELECTOR_KINDS.md decision #16. The
    # lowering must produce a 2D `index` op whose operand is `\$Fx` for
    # face_x entries and `\$Fy` for face_y entries; axis-args are the sorted
    # axis pattern variables (`\$x` first, then `\$y`), with the entry's
    # axis slot carrying the offset. Stagger → operand mapping: sorted unique
    # staggers (face_x < face_y) correspond to args in order ($Fx, $Fy).
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "divergence_arakawa_c.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["divergence_arakawa_c"])
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 4

    # Entry 1: face_x, axis $x, offset 0 -> index($Fx, $x, $y), coeff -1/dx
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test length(idx1["args"]) == 3
    @test String(idx1["args"][1]) == "\$Fx"
    @test String(idx1["args"][2]) == "\$x"   # offset == 0 -> bare axis
    @test String(idx1["args"][3]) == "\$y"

    # Entry 2: face_x, axis $x, offset 1 -> index($Fx, $x + 1, $y)
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test length(idx2["args"]) == 3
    @test String(idx2["args"][1]) == "\$Fx"
    a2_x = idx2["args"][2]
    @test a2_x["op"] == "+"
    @test String(a2_x["args"][1]) == "\$x"
    @test Int(a2_x["args"][2]) == 1
    @test String(idx2["args"][3]) == "\$y"

    # Entry 3: face_y, axis $y, offset 0 -> index($Fy, $x, $y)
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    @test length(idx3["args"]) == 3
    @test String(idx3["args"][1]) == "\$Fy"
    @test String(idx3["args"][2]) == "\$x"
    @test String(idx3["args"][3]) == "\$y"

    # Entry 4: face_y, axis $y, offset 1 -> index($Fy, $x, $y + 1)
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    @test length(idx4["args"]) == 3
    @test String(idx4["args"][1]) == "\$Fy"
    @test String(idx4["args"][2]) == "\$x"
    a4_y = idx4["args"][3]
    @test a4_y["op"] == "+"
    @test String(a4_y["args"][1]) == "\$y"
    @test Int(a4_y["args"][2]) == 1

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

@testitem "lower_stencil_to_replacement: arakawa errors on bad stagger" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[Dict(
            "selector" => Dict(
                "kind" => "arakawa",
                "stagger" => "wat",
                "axis" => "\$x",
                "offset" => 0,
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
    @test occursin("'wat'", err.msg)
end

@testitem "lower_stencil_to_replacement: arakawa errors on non-pattern-var axis" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[Dict(
            "selector" => Dict(
                "kind" => "arakawa",
                "stagger" => "face_x",
                "axis" => "x",
                "offset" => 0,
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
    @test occursin("\$", err.msg)
end

@testitem "lower_stencil_to_replacement: arakawa works without applies_to.dim" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement

    # Arakawa rules carry no top-level `applies_to.dim` — axis pattern
    # variables are intrinsic to each stencil entry. The lowerer must
    # accept this shape (which the cartesian path requires `dim` for).
    rule = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[
            Dict(
                "selector" => Dict(
                    "kind" => "arakawa", "stagger" => "face_x",
                    "axis" => "\$x", "offset" => 0,
                ),
                "coeff" => -1,
            ),
            Dict(
                "selector" => Dict(
                    "kind" => "arakawa", "stagger" => "face_x",
                    "axis" => "\$x", "offset" => 1,
                ),
                "coeff" => 1,
            ),
        ],
    )
    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")
    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2
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

@testitem "lower_stencil_to_replacement: lowers lax_friedrichs_flux (cartesian)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_volume",
        "lax_friedrichs_flux.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["lax_friedrichs_flux"])
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")
    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # First operand of `applies_to.args = ["\$q", "\$c"]` is `\$q` —
    # the stencil indexes q, not c.
    e1 = repl["args"][1]
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$q"
end

@testitem "lower_stencil_to_replacement: lowers centered_2nd_deriv_uniform (cartesian)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "centered_2nd_deriv_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["centered_2nd_deriv_uniform"])
    @test !haskey(rule, "replacement")
    @test rule["applies_to"]["op"] == "d2"
    @test String(rule["applies_to"]["args"][1]) == "\$u"
    @test String(rule["applies_to"]["dim"]) == "\$x"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    # Three stencil entries → combine(term1, term2, term3)
    @test repl["op"] == "+"
    @test length(repl["args"]) == 3

    # Entry 1: offset -1 → coeff * index($u, $x + (-1))
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    arg1 = idx1["args"][2]
    @test arg1["op"] == "+"
    @test String(arg1["args"][1]) == "\$x"
    @test Int(arg1["args"][2]) == -1

    # Entry 2: offset 0 → coeff * index($u, $x) — no `+ 0` wrapper
    e2 = repl["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$x"

    # Entry 3: offset +1 → coeff * index($u, $x + 1)
    e3 = repl["args"][3]
    @test e3["op"] == "*"
    idx3 = e3["args"][2]
    @test idx3["op"] == "index"
    @test String(idx3["args"][1]) == "\$u"
    arg3 = idx3["args"][2]
    @test arg3["op"] == "+"
    @test String(arg3["args"][1]) == "\$x"
    @test Int(arg3["args"][2]) == 1

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing

    # Idempotence: lowering again returns the same replacement
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: lowers laplacian_2nd_uniform_cartesian (arakawa)" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "laplacian_2nd_uniform_cartesian.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["laplacian_2nd_uniform_cartesian"])
    @test !haskey(rule, "replacement")
    @test rule["applies_to"]["op"] == "laplacian"
    @test String(rule["applies_to"]["args"][1]) == "\$u"

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    # 6 stencil entries (3 for $x, 3 for $y) → combine(6 terms)
    @test repl["op"] == "+"
    @test length(repl["args"]) == 6

    # All terms are coeff * index($u, $x_arg, $y_arg)
    for (i, term) in enumerate(repl["args"])
        @test term["op"] == "*"
        idx = term["args"][2]
        @test idx["op"] == "index"
        @test length(idx["args"]) == 3   # operand + two axis args ($x, $y)
        @test String(idx["args"][1]) == "\$u"
    end

    # Entries 1–3 shift $x: each index has $x+offset or bare $x in slot 2, bare $y in slot 3
    # Entry 1: axis=$x offset=-1 → index($u, $x+(-1), $y)
    e1 = repl["args"][1]
    idx1 = e1["args"][2]
    @test idx1["args"][2]["op"] == "+"       # $x + (-1)
    @test String(idx1["args"][2]["args"][1]) == "\$x"
    @test Int(idx1["args"][2]["args"][2]) == -1
    @test String(idx1["args"][3]) == "\$y"   # $y bare

    # Entry 2: axis=$x offset=0 → index($u, $x, $y) — no `+ 0` wrapper
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test String(idx2["args"][2]) == "\$x"
    @test String(idx2["args"][3]) == "\$y"

    # Entry 3: axis=$x offset=+1 → index($u, $x+1, $y)
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    @test idx3["args"][2]["op"] == "+"
    @test Int(idx3["args"][2]["args"][2]) == 1
    @test String(idx3["args"][3]) == "\$y"

    # Entries 4–6 shift $y: each index has bare $x in slot 2, $y+offset or bare $y in slot 3
    # Entry 4: axis=$y offset=-1 → index($u, $x, $y+(-1))
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    @test String(idx4["args"][2]) == "\$x"   # $x bare
    @test idx4["args"][3]["op"] == "+"        # $y + (-1)
    @test String(idx4["args"][3]["args"][1]) == "\$y"
    @test Int(idx4["args"][3]["args"][2]) == -1

    # Entry 5: axis=$y offset=0 → index($u, $x, $y) — no `+ 0` wrapper
    e5 = repl["args"][5]
    idx5 = e5["args"][2]
    @test String(idx5["args"][2]) == "\$x"
    @test String(idx5["args"][3]) == "\$y"

    # Entry 6: axis=$y offset=+1 → index($u, $x, $y+1)
    e6 = repl["args"][6]
    idx6 = e6["args"][2]
    @test String(idx6["args"][2]) == "\$x"
    @test idx6["args"][3]["op"] == "+"
    @test Int(idx6["args"][3]["args"][2]) == 1

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(out["replacement"])
    @test expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: lowers staggered_1st_uniform_cc_to_face" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # staggered_1st_uniform_cc_to_face: arakawa rule with cell_center selectors,
    # 1D x-axis, offsets {-1, 0}. Mirrors MOL StaggeredStencilInfo
    # EdgeAlignedVar: d/dx(u)|_face[j] = (u_cc[j] - u_cc[j-1]) / h.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "staggered_1st_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["staggered_1st_uniform_cc_to_face"])
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # Entry 1: cell_center, axis $x, offset -1 -> index($u, $x + (-1)), coeff -1/h
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    arg1 = idx1["args"][2]
    @test arg1["op"] == "+"
    @test String(arg1["args"][1]) == "\$x"
    @test Int(arg1["args"][2]) == -1

    # Entry 2: cell_center, axis $x, offset 0 -> index($u, $x) (no `+ 0` wrapper)
    e2 = repl["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    @test String(idx2["args"][2]) == "\$x"

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(repl)
    @test expr !== nothing

    # Idempotence
    @test lower_stencil_to_replacement(out)["replacement"] == repl
end

@testitem "lower_stencil_to_replacement: lowers staggered_1st_uniform_face_to_cc" begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    import EarthSciSerialization
    using JSON

    # staggered_1st_uniform_face_to_cc: arakawa rule with face_x selectors,
    # 1D x-axis, offsets {0, +1}. Mirrors MOL StaggeredStencilInfo
    # CenterAlignedVar: d/dx(u)|_cc[i] = (u_face[i+1] - u_face[i]) / h.
    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "staggered_1st_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["staggered_1st_uniform_face_to_cc"])
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
    @test haskey(out, "replacement")

    repl = out["replacement"]
    @test repl["op"] == "+"
    @test length(repl["args"]) == 2

    # Entry 1: face_x, axis $x, offset 0 -> index($u, $x) (no `+ 0` wrapper), coeff -1/h
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    @test String(idx1["args"][2]) == "\$x"

    # Entry 2: face_x, axis $x, offset 1 -> index($u, $x + 1), coeff +1/h
    e2 = repl["args"][2]
    @test e2["op"] == "*"
    idx2 = e2["args"][2]
    @test idx2["op"] == "index"
    @test String(idx2["args"][1]) == "\$u"
    arg2 = idx2["args"][2]
    @test arg2["op"] == "+"
    @test String(arg2["args"][1]) == "\$x"
    @test Int(arg2["args"][2]) == 1

    # ESS parse_expression accepts the lowered AST
    expr = EarthSciSerialization.parse_expression(repl)
    @test expr !== nothing

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

@testitem "lower_stencil_to_scheme: lowers cartesian upwind_1st to ESS scheme parts" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st"])

    scheme, use_rule = lower_stencil_to_scheme("upwind_1st", rule)

    # Whitelisted structural subset only — the catalog's sibling
    # `replacement` AST must not leak into the ESS scheme object.
    @test Set(keys(scheme)) == Set(["applies_to", "grid_family", "combine", "stencil", "accuracy"])
    @test scheme["grid_family"] == "cartesian"
    @test scheme["combine"] == "+"
    @test scheme["applies_to"] == rule["applies_to"]
    @test length(scheme["stencil"]) == 2
    sel1 = scheme["stencil"][1]["selector"]
    @test sel1 == Dict{String, Any}("kind" => "cartesian", "axis" => "\$x", "offset" => -1)
    @test scheme["stencil"][1]["coeff"] == rule["stencil"][1]["coeff"]

    # The use: rule re-states applies_to as its pattern so pattern variables
    # bind one-to-one to the scheme operands.
    @test use_rule == Dict{String, Any}(
        "name"    => "upwind_1st",
        "pattern" => rule["applies_to"],
        "use"     => "upwind_1st",
    )
end

@testitem "lower_stencil_to_scheme: ESS parse_schemes + parse_rule accept the lowered parts" begin
    using EarthSciDiscretizations: lower_stencil_to_scheme
    import EarthSciSerialization
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "upwind_1st.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["upwind_1st"])

    scheme, use_rule = lower_stencil_to_scheme("upwind_1st", rule)
    schemes = EarthSciSerialization.parse_schemes(Dict{String, Any}("upwind_1st" => scheme))
    @test haskey(schemes, "upwind_1st")
    parsed_rule = EarthSciSerialization.parse_rule(use_rule)
    @test parsed_rule.replacement_scheme == "upwind_1st"
end

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
    @test occursin("cartesian-only", err.msg)

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

@testitem "lower_stencil_to_canonical_replacement: multi-axis laplacian → i/j components" begin
    using EarthSciDiscretizations: lower_stencil_to_canonical_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "laplacian_2nd_uniform_cartesian.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["laplacian_2nd_uniform_cartesian"])

    expr = lower_stencil_to_canonical_replacement(rule)
    @test expr["op"] == "+"
    @test length(expr["args"]) == 6

    # No axis pattern variables survive; the operand pattern variable does.
    s = JSON.json(expr)
    @test !occursin("\$x", s)
    @test !occursin("\$y", s)
    @test occursin("\$u", s)

    # Entry 1: $x-axis offset -1 with $y untouched → index($u, i - 1, j)
    # (sorted axis order: $x → i, $y → j; arakawa lowering keeps the
    # entry's own axis offset and the bare component elsewhere).
    idx1 = expr["args"][1]["args"][2]
    @test idx1["op"] == "index"
    @test String(idx1["args"][1]) == "\$u"
    @test idx1["args"][2]["op"] == "+"
    @test String(idx1["args"][2]["args"][1]) == "i"
    @test Int(idx1["args"][2]["args"][2]) == -1
    @test String(idx1["args"][3]) == "j"
end

@testitem "lower_stencil_to_canonical_replacement: replacement-form rule → dim var becomes i" begin
    using EarthSciDiscretizations: lower_stencil_to_canonical_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "centered_2nd_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["centered_2nd_uniform"])

    expr = lower_stencil_to_canonical_replacement(rule)
    s = JSON.json(expr)
    # The arrayop wrapper is stripped and $x is substituted with i.
    @test expr["op"] != "arrayop"
    @test !occursin("\$x", s)
    @test occursin("\$u", s)
    @test occursin("\"i\"", s)
end

@testitem "lower_stencil_to_canonical_replacement: rejects staggered and literal-axis stencils" begin
    using EarthSciDiscretizations: lower_stencil_to_canonical_replacement

    staggered = Dict{String, Any}(
        "applies_to" => Dict("op" => "div", "args" => ["\$F"]),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "arakawa", "axis" => "\$x", "offset" => 0,
                               "stagger" => "face_x"),
            "coeff" => 1,
        )],
    )
    err = try lower_stencil_to_canonical_replacement(staggered); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("cell_center", err.msg)

    literal_axis = Dict{String, Any}(
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "lat"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "latlon", "axis" => "lat", "offset" => 1),
            "coeff" => 1,
        )],
    )
    err = try lower_stencil_to_canonical_replacement(literal_axis); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("pattern variable", err.msg)
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

@testitem "lower_stencil_to_canonical_replacement: nested-dim pattern maps every axis" begin
    using EarthSciDiscretizations: lower_stencil_to_canonical_replacement
    using JSON

    path = joinpath(
        dirname(dirname(@__FILE__)),
        "discretizations",
        "finite_difference",
        "mixed_deriv_2nd_uniform.json",
    )
    raw = JSON.parsefile(path)
    rule = Dict{String, Any}(raw["discretizations"]["mixed_deriv_2nd_uniform"])

    # The pattern is grad(grad($u, dim=$y), dim=$x): both nested dim pattern
    # variables must map to canonical components (sorted: $x → i, $y → j).
    expr = lower_stencil_to_canonical_replacement(rule)
    s = JSON.json(expr)
    @test !occursin("\$x", s)
    @test !occursin("\$y", s)
    @test occursin("\$u", s)
    @test occursin("\"i\"", s)
    @test occursin("\"j\"", s)
end
