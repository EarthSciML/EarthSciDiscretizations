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
    @test !haskey(rule, "replacement")

    out = lower_stencil_to_replacement(rule)
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
        "applies_to" => Dict("op" => "grad", "args" => ["\$u"], "dim" => "\$k"),
        "stencil" => Any[Dict(
            "selector" => Dict("kind" => "latlon", "axis" => "lon", "offset" => 1),
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
    @test occursin("'latlon'", err.msg)
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

    # `divergence_arakawa_c` is the canonical arakawa-kind rule today
    # (per SELECTOR_KINDS.md decision #13 / #16). The lowering must
    # produce a 2D `index` op on the operand `\$F` whose positional
    # axis-args are the sorted axis pattern variables (`\$x` first,
    # then `\$y`), with the entry's axis slot carrying the offset.
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

    # Entry 1: face_x, axis $x, offset 0 -> index($F, $x, $y), coeff -1/dx
    e1 = repl["args"][1]
    @test e1["op"] == "*"
    idx1 = e1["args"][2]
    @test idx1["op"] == "index"
    @test length(idx1["args"]) == 3
    @test String(idx1["args"][1]) == "\$F"
    @test String(idx1["args"][2]) == "\$x"   # offset == 0 -> bare axis
    @test String(idx1["args"][3]) == "\$y"

    # Entry 2: face_x, axis $x, offset 1 -> index($F, $x + 1, $y)
    e2 = repl["args"][2]
    idx2 = e2["args"][2]
    @test length(idx2["args"]) == 3
    a2_x = idx2["args"][2]
    @test a2_x["op"] == "+"
    @test String(a2_x["args"][1]) == "\$x"
    @test Int(a2_x["args"][2]) == 1
    @test String(idx2["args"][3]) == "\$y"

    # Entry 3: face_y, axis $y, offset 0 -> index($F, $x, $y)
    e3 = repl["args"][3]
    idx3 = e3["args"][2]
    @test length(idx3["args"]) == 3
    @test String(idx3["args"][2]) == "\$x"
    @test String(idx3["args"][3]) == "\$y"

    # Entry 4: face_y, axis $y, offset 1 -> index($F, $x, $y + 1)
    e4 = repl["args"][4]
    idx4 = e4["args"][2]
    @test length(idx4["args"]) == 3
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
