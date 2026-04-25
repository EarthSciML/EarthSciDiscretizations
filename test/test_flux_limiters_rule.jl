using Test
using TestItems

# Tests for the finite_volume/flux_limiter_{minmod,superbee} declarative rules.
#
# Both rules encode the limiter formula phi(r) as a fully AST-encoded
# ExpressionNode (ESS §7 §9.2) — there is no runtime-callable Julia name.
# The tests exercise:
#
# Layer A   — rule discovery + JSON byte-diff round-trip.
# Layer A'  — AST op validation: every op in the formula is one of the ESS
#             ExpressionNode ops and the tree is composed of Dict/String/Number
#             nodes only.
# Layer B   — AST-evaluated phi(r) via a minimal pure-Julia ExpressionNode
#             evaluator scoped to this test harness. ESS 0.0.3's evaluator
#             does not yet implement the `max`/`min` ops required by these
#             limiters (see follow-up bead filed against ESS); the bead risk
#             section explicitly allows a local evaluator here. Verifies
#             positivity, monotonicity-preserving behavior (phi(r)=0 for
#             r<=0), Sweby upper bound phi(r)<=2, consistency phi(1)=1, and
#             exact agreement with fixture reference (r, phi) pairs.
# Layer B'  — 1D periodic linear advection with a MUSCL-style upwind scheme
#             whose slope correction is scaled by phi(r) evaluated from the
#             rule's AST on every interface at every step. Verifies the
#             strict TVD property TV(q_final) <= TV(q_initial) + tol.

# ---------------------------------------------------------------------------
# Layer A: discovery + byte-diff
# ---------------------------------------------------------------------------

@testitem "flux_limiter_minmod rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "flux_limiter_minmod", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"limit\"", content)
    @test occursin("\"tvd_flux_limiter\"", content)
    @test occursin("\"formula\"", content)
    @test occursin("\"op\": \"max\"", content)
    @test occursin("\"op\": \"min\"", content)
    # The minmod formula must contain the slope-ratio variable placeholder.
    @test occursin("\"\$r\"", content)
end

@testitem "flux_limiter_superbee rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "flux_limiter_superbee", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"limit\"", content)
    @test occursin("\"tvd_flux_limiter\"", content)
    @test occursin("\"formula\"", content)
    @test occursin("\"op\": \"max\"", content)
    @test occursin("\"op\": \"min\"", content)
    @test occursin("\"op\": \"*\"", content)
    @test occursin("\"\$r\"", content)
end

@testitem "flux limiter rule JSON round-trips byte-stable (both limiters)" begin
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    for name in ("flux_limiter_minmod", "flux_limiter_superbee")
        idx = findfirst(r -> r.name == name, rules)
        @test idx !== nothing
        rule = rules[idx]

        raw = read(rule.path, String)
        parsed = JSON.parse(raw)
        @test parsed isa AbstractDict
        @test haskey(parsed, "discretizations")
        @test haskey(parsed["discretizations"], name)

        reserialized = JSON.json(parsed)
        reparsed = JSON.parse(reserialized)
        @test reparsed == parsed

        spec = parsed["discretizations"][name]
        @test spec["applies_to"]["op"] == "limit"
        @test spec["applies_to"]["args"] == ["\$r"]
        @test spec["grid_family"] == "cartesian"
        @test spec["form"] == "tvd_flux_limiter"
        @test spec["tvd"] == true
        @test spec["monotonicity_preserving"] == true
        @test haskey(spec, "formula")
        @test haskey(spec["formula"], "op")
        @test haskey(spec["formula"], "args")
    end
end

# ---------------------------------------------------------------------------
# Layer A': AST op validation
# ---------------------------------------------------------------------------

@testitem "flux limiter AST uses only ESS ExpressionNode ops" begin
    using EarthSciDiscretizations: load_rules
    using JSON

    # Ops declared by ESS §7 ExpressionNode $def (esm-spec.md line 126, 134).
    allowed_ops = Set(
        [
            "+", "-", "*", "/", "^",
            "max", "min", "ifelse", "abs", "sign",
            "exp", "log", "log10", "sqrt",
            "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
            "floor", "ceil",
            "==", "!=", "<", "<=", ">", ">=",
            "and", "or", "not",
        ]
    )

    function collect_ops!(acc, node)
        if node isa AbstractDict
            if haskey(node, "op")
                push!(acc, node["op"])
            end
            if haskey(node, "args")
                for a in node["args"]
                    collect_ops!(acc, a)
                end
            end
        elseif node isa AbstractVector
            for a in node
                collect_ops!(acc, a)
            end
        end
        return acc
    end

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    for name in ("flux_limiter_minmod", "flux_limiter_superbee")
        idx = findfirst(r -> r.name == name, rules)
        rule = rules[idx]
        parsed = JSON.parse(read(rule.path, String))
        formula = parsed["discretizations"][name]["formula"]
        ops = collect_ops!(String[], formula)
        @test !isempty(ops)
        for op in ops
            @test op in allowed_ops
        end
    end
end

# ---------------------------------------------------------------------------
# Layer B: AST-evaluated phi(r) against fixture reference pairs
# ---------------------------------------------------------------------------

@testitem "flux_limiter_minmod AST reproduces reference phi(r) values" begin
    using EarthSciDiscretizations
    using JSON

    # Minimal pure-Julia evaluator scoped to the test harness. ESS 0.0.3's
    # `evaluate` does not yet implement the `max`/`min` ops used by these
    # limiters (see bead risk section). Supports: Number literals, String
    # variables of form "$name" (resolved against bindings by stripping the
    # leading "$"), and Dict nodes with ops max, min, +, -, *, /.
    function eval_ast(node, bindings::Dict{String, Float64})::Float64
        if node isa Number
            return Float64(node)
        elseif node isa AbstractString
            name = startswith(node, "\$") ? node[2:end] : node
            haskey(bindings, name) ||
                throw(ArgumentError("unbound variable: $name"))
            return bindings[name]
        elseif node isa AbstractDict
            op = node["op"]
            args = [eval_ast(a, bindings) for a in node["args"]]
            if op == "max"
                return maximum(args)
            elseif op == "min"
                return minimum(args)
            elseif op == "+"
                return length(args) == 1 ? args[1] : sum(args)
            elseif op == "-"
                return length(args) == 1 ? -args[1] : args[1] - args[2]
            elseif op == "*"
                return length(args) == 1 ? args[1] : prod(args)
            elseif op == "/"
                return args[1] / args[2]
            else
                throw(ArgumentError("unsupported op in test-harness eval: $op"))
            end
        end
        throw(ArgumentError("unrecognized AST node type: $(typeof(node))"))
    end

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_volume", "flux_limiter_minmod.json")
    fixture = JSON.parse(
        read(
            joinpath(repo_root, "tests", "fixtures", "flux_limiter_minmod", "phi_sweep.esm"),
            String,
        )
    )
    @test fixture["rule"] == "flux_limiter_minmod"

    spec = JSON.parse(read(rule_path, String))["discretizations"]["flux_limiter_minmod"]
    formula = spec["formula"]
    tol = Float64(fixture["tolerance"])

    for pair in fixture["reference_values"]
        r = Float64(pair["r"])
        phi_expected = Float64(pair["phi"])
        phi_actual = eval_ast(formula, Dict("r" => r))
        @test isapprox(phi_actual, phi_expected; atol = tol, rtol = 0)
    end

    sweep_min = Float64(fixture["tvd_properties"]["sweep_r_min"])
    sweep_max = Float64(fixture["tvd_properties"]["sweep_r_max"])
    sweep_step = Float64(fixture["tvd_properties"]["sweep_r_step"])
    for r in sweep_min:sweep_step:sweep_max
        phi = eval_ast(formula, Dict("r" => Float64(r)))
        @test phi >= -tol
        @test phi <= 2.0 + tol
        if r <= 0
            @test isapprox(phi, 0.0; atol = tol)
        end
    end

    @test isapprox(eval_ast(formula, Dict("r" => 1.0)), 1.0; atol = tol)
end

@testitem "flux_limiter_superbee AST reproduces reference phi(r) values" begin
    using EarthSciDiscretizations
    using JSON

    function eval_ast(node, bindings::Dict{String, Float64})::Float64
        if node isa Number
            return Float64(node)
        elseif node isa AbstractString
            name = startswith(node, "\$") ? node[2:end] : node
            haskey(bindings, name) ||
                throw(ArgumentError("unbound variable: $name"))
            return bindings[name]
        elseif node isa AbstractDict
            op = node["op"]
            args = [eval_ast(a, bindings) for a in node["args"]]
            if op == "max"
                return maximum(args)
            elseif op == "min"
                return minimum(args)
            elseif op == "+"
                return length(args) == 1 ? args[1] : sum(args)
            elseif op == "-"
                return length(args) == 1 ? -args[1] : args[1] - args[2]
            elseif op == "*"
                return length(args) == 1 ? args[1] : prod(args)
            elseif op == "/"
                return args[1] / args[2]
            else
                throw(ArgumentError("unsupported op in test-harness eval: $op"))
            end
        end
        throw(ArgumentError("unrecognized AST node type: $(typeof(node))"))
    end

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_volume", "flux_limiter_superbee.json")
    fixture = JSON.parse(
        read(
            joinpath(repo_root, "tests", "fixtures", "flux_limiter_superbee", "phi_sweep.esm"),
            String,
        )
    )
    @test fixture["rule"] == "flux_limiter_superbee"

    spec = JSON.parse(read(rule_path, String))["discretizations"]["flux_limiter_superbee"]
    formula = spec["formula"]
    tol = Float64(fixture["tolerance"])

    for pair in fixture["reference_values"]
        r = Float64(pair["r"])
        phi_expected = Float64(pair["phi"])
        phi_actual = eval_ast(formula, Dict("r" => r))
        @test isapprox(phi_actual, phi_expected; atol = tol, rtol = 0)
    end

    sweep_min = Float64(fixture["tvd_properties"]["sweep_r_min"])
    sweep_max = Float64(fixture["tvd_properties"]["sweep_r_max"])
    sweep_step = Float64(fixture["tvd_properties"]["sweep_r_step"])
    for r in sweep_min:sweep_step:sweep_max
        phi = eval_ast(formula, Dict("r" => Float64(r)))
        @test phi >= -tol
        @test phi <= 2.0 + tol
        if r <= 0
            @test isapprox(phi, 0.0; atol = tol)
        end
    end

    @test isapprox(eval_ast(formula, Dict("r" => 1.0)), 1.0; atol = tol)
end

# ---------------------------------------------------------------------------
# Layer B': 1D periodic advection TVD check via AST-evaluated limiter
# ---------------------------------------------------------------------------

@testitem "flux_limiter_minmod: 1D advection TVD property holds" begin
    using EarthSciDiscretizations
    using JSON

    function eval_ast(node, bindings::Dict{String, Float64})::Float64
        if node isa Number
            return Float64(node)
        elseif node isa AbstractString
            name = startswith(node, "\$") ? node[2:end] : node
            haskey(bindings, name) ||
                throw(ArgumentError("unbound variable: $name"))
            return bindings[name]
        elseif node isa AbstractDict
            op = node["op"]
            args = [eval_ast(a, bindings) for a in node["args"]]
            if op == "max"
                return maximum(args)
            elseif op == "min"
                return minimum(args)
            elseif op == "+"
                return length(args) == 1 ? args[1] : sum(args)
            elseif op == "-"
                return length(args) == 1 ? -args[1] : args[1] - args[2]
            elseif op == "*"
                return length(args) == 1 ? args[1] : prod(args)
            elseif op == "/"
                return args[1] / args[2]
            else
                throw(ArgumentError("unsupported op: $op"))
            end
        end
        throw(ArgumentError("unrecognized AST node type: $(typeof(node))"))
    end

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_volume", "flux_limiter_minmod.json")
    fixture = JSON.parse(
        read(
            joinpath(repo_root, "tests", "fixtures", "flux_limiter_minmod", "tvd_check.esm"),
            String,
        )
    )

    spec = JSON.parse(read(rule_path, String))["discretizations"]["flux_limiter_minmod"]
    formula = spec["formula"]

    n = Int(fixture["grid"]["n"])
    dx = Float64(fixture["grid"]["dx"])
    u = Float64(fixture["advection"]["velocity"])
    cfl = Float64(fixture["advection"]["cfl"])
    periods = Float64(fixture["advection"]["periods"])
    eps_denom = Float64(fixture["eps_denom"])
    tvd_tol = Float64(fixture["tvd_tolerance"])

    # Analytical cell averages of the smooth+square initial condition. Smooth
    # component: sin(2*pi*x) on [0, 0.4]; square component: 1 on [0.6, 0.85].
    F_smooth(x) = -cos(2pi * x) / (2pi)
    function cell_avg(i)
        a = (i - 1) * dx
        b = i * dx
        sa = max(a, 0.0); sb = min(b, 0.4)
        smooth = sb > sa ? (F_smooth(sb) - F_smooth(sa)) / dx : 0.0
        qa = max(a, 0.6); qb = min(b, 0.85)
        square = qb > qa ? (qb - qa) / dx : 0.0
        return smooth + square
    end
    q0 = [cell_avg(i) for i in 1:n]

    modn(j) = mod(j - 1, n) + 1
    tv(q) = sum(abs(q[modn(i + 1)] - q[i]) for i in 1:n)

    phi(r) = eval_ast(formula, Dict("r" => Float64(r)))

    function step!(q, dt)
        Fs = Vector{Float64}(undef, n)
        for i in 1:n
            qm = q[modn(i - 1)]; qi = q[i]; qp = q[modn(i + 1)]
            num = qi - qm; den = qp - qi
            den_safe = abs(den) > eps_denom ? den : (den >= 0 ? eps_denom : -eps_denom)
            r_ratio = num / den_safe
            Fs[i] = u * (qi + 0.5 * phi(r_ratio) * (qp - qi))
        end
        qnew = similar(q)
        for i in 1:n
            qnew[i] = q[i] - dt / dx * (Fs[i] - Fs[modn(i - 1)])
        end
        return qnew
    end

    # Wrap the time loop in a function so loop-local variables live in hard
    # scope (testitem module-scope `while` otherwise hits soft-scope ambiguity).
    function advance(q0_in, dt, T)
        q = copy(q0_in)
        tsim = 0.0
        while tsim < T
            dts = min(dt, T - tsim)
            q = step!(q, dts)
            tsim += dts
        end
        return q
    end

    dt = cfl * dx / abs(u)
    T = periods / abs(u)
    q = advance(q0, dt, T)

    tv_initial = tv(q0)
    tv_final = tv(q)
    @info "flux_limiter_minmod TVD" tv_initial tv_final
    @test tv_final <= tv_initial + tvd_tol
    @test all(isfinite, q)
end

@testitem "flux_limiter_superbee: 1D advection TVD property holds" begin
    using EarthSciDiscretizations
    using JSON

    function eval_ast(node, bindings::Dict{String, Float64})::Float64
        if node isa Number
            return Float64(node)
        elseif node isa AbstractString
            name = startswith(node, "\$") ? node[2:end] : node
            haskey(bindings, name) ||
                throw(ArgumentError("unbound variable: $name"))
            return bindings[name]
        elseif node isa AbstractDict
            op = node["op"]
            args = [eval_ast(a, bindings) for a in node["args"]]
            if op == "max"
                return maximum(args)
            elseif op == "min"
                return minimum(args)
            elseif op == "+"
                return length(args) == 1 ? args[1] : sum(args)
            elseif op == "-"
                return length(args) == 1 ? -args[1] : args[1] - args[2]
            elseif op == "*"
                return length(args) == 1 ? args[1] : prod(args)
            elseif op == "/"
                return args[1] / args[2]
            else
                throw(ArgumentError("unsupported op: $op"))
            end
        end
        throw(ArgumentError("unrecognized AST node type: $(typeof(node))"))
    end

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_volume", "flux_limiter_superbee.json")
    fixture = JSON.parse(
        read(
            joinpath(repo_root, "tests", "fixtures", "flux_limiter_superbee", "tvd_check.esm"),
            String,
        )
    )

    spec = JSON.parse(read(rule_path, String))["discretizations"]["flux_limiter_superbee"]
    formula = spec["formula"]

    n = Int(fixture["grid"]["n"])
    dx = Float64(fixture["grid"]["dx"])
    u = Float64(fixture["advection"]["velocity"])
    cfl = Float64(fixture["advection"]["cfl"])
    periods = Float64(fixture["advection"]["periods"])
    eps_denom = Float64(fixture["eps_denom"])
    tvd_tol = Float64(fixture["tvd_tolerance"])

    F_smooth(x) = -cos(2pi * x) / (2pi)
    function cell_avg(i)
        a = (i - 1) * dx
        b = i * dx
        sa = max(a, 0.0); sb = min(b, 0.4)
        smooth = sb > sa ? (F_smooth(sb) - F_smooth(sa)) / dx : 0.0
        qa = max(a, 0.6); qb = min(b, 0.85)
        square = qb > qa ? (qb - qa) / dx : 0.0
        return smooth + square
    end
    q0 = [cell_avg(i) for i in 1:n]

    modn(j) = mod(j - 1, n) + 1
    tv(q) = sum(abs(q[modn(i + 1)] - q[i]) for i in 1:n)

    phi(r) = eval_ast(formula, Dict("r" => Float64(r)))

    function step!(q, dt)
        Fs = Vector{Float64}(undef, n)
        for i in 1:n
            qm = q[modn(i - 1)]; qi = q[i]; qp = q[modn(i + 1)]
            num = qi - qm; den = qp - qi
            den_safe = abs(den) > eps_denom ? den : (den >= 0 ? eps_denom : -eps_denom)
            r_ratio = num / den_safe
            Fs[i] = u * (qi + 0.5 * phi(r_ratio) * (qp - qi))
        end
        qnew = similar(q)
        for i in 1:n
            qnew[i] = q[i] - dt / dx * (Fs[i] - Fs[modn(i - 1)])
        end
        return qnew
    end

    # Wrap the time loop in a function so loop-local variables live in hard
    # scope (testitem module-scope `while` otherwise hits soft-scope ambiguity).
    function advance(q0_in, dt, T)
        q = copy(q0_in)
        tsim = 0.0
        while tsim < T
            dts = min(dt, T - tsim)
            q = step!(q, dts)
            tsim += dts
        end
        return q
    end

    dt = cfl * dx / abs(u)
    T = periods / abs(u)
    q = advance(q0, dt, T)

    tv_initial = tv(q0)
    tv_final = tv(q)
    @info "flux_limiter_superbee TVD" tv_initial tv_final
    @test tv_final <= tv_initial + tvd_tol
    @test all(isfinite, q)
end

# ---------------------------------------------------------------------------
# Catalog-level smoke: both limiter rules appear under :finite_volume.
# ---------------------------------------------------------------------------

@testitem "flux limiter rules are present in the :finite_volume family" begin
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    fv = Set(r.name for r in rules if r.family == :finite_volume)
    @test "flux_limiter_minmod" in fv
    @test "flux_limiter_superbee" in fv
end
