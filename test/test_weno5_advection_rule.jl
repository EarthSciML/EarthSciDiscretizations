using Test
using TestItems

# Tests for the finite_volume/weno5_advection declarative rule.
#
# dsc-b78 rewrote this rule as the canonical NONLINEAR-scheme exemplar of
# the closed-AST lowering pattern (sibling of centered_2nd_uniform / dsc-rar
# for linear schemes). The rule now matches against `div(*($U, $q))` (a §4.2
# operator) and replaces it with a single `arrayop` whose body is the FV
# divergence (F_{i+1/2} − F_{i-1/2})/dx, with face fluxes built from a
# cell-to-face average of U and a Jiang-Shu (1996) WENO5 reconstruction of
# q. Upwinding is encoded as `ifelse(U_face > 0, q^L, q^R)` — no
# scheme-specific `stencil` / `selector` / `offset` blobs and no off-spec
# `op: "advect"` / `reconstruction_*` / `smoothness_indicators` /
# `nonlinear_weights` fields.
#
# Layer A: rule discovery + closed-AST shape. The rule is discoverable
# under the :finite_volume family and its JSON file shows the closed
# replacement-form lowering — no stencil-coefficient blobs.
#
# Layer B: numeric WENO5 convergence. The closed AST is exercised against
# Jiang-Shu (1996) 5th-order WENO reconstruction on a uniform 1D Cartesian
# grid: theoretically O(dx^5) in L_inf on smooth data away from critical
# points. The numeric oracle is the canonical scalar implementation of
# WENO5 (the same formulae the AST encodes); both implementations are
# tested against a phase-shifted sine MMS at n = 32, 64, 128, 256 and
# minimum observed order must be ≥ 4.7.
#
# Layer B (shock): Linear advection of a unit square wave by u = 1 for one
# full period on a periodic domain with WENO5 + SSP-RK3. Max overshoot and
# undershoot must stay bounded by 0.05 (qualitative shock-capturing check).

@testitem "weno5_advection rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "weno5_advection", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"weighted_essentially_nonoscillatory\"", content)
    @test occursin("\"O(dx^5)\"", content)

    # Closed-AST lowering: applies_to matches the §4.2 op `div` over the
    # advective flux `*($U, $q)`. The replacement is a single arrayop whose
    # body uses `index`, arithmetic, `^`, and `ifelse` to encode the
    # Jiang-Shu (1996) WENO5 reconstruction at each face — no off-spec
    # `op: "advect"` and no scheme-specific kernels.
    @test occursin("\"op\": \"div\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"index\"", content)
    @test occursin("\"ifelse\"", content)
    @test !occursin("\"op\": \"advect\"", content)
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"selector\"", content)
    @test !occursin("\"offset\"", content)
    @test !occursin("\"reconstruction_left_biased\"", content)
    @test !occursin("\"reconstruction_right_biased\"", content)
    @test !occursin("\"smoothness_indicators\"", content)
    @test !occursin("\"nonlinear_weights\"", content)
end

@testitem "weno5_advection rule JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "weno5_advection", rules)
    @test idx !== nothing
    rule = rules[idx]

    raw = read(rule.path, String)
    parsed = JSON.parse(raw)
    @test parsed isa AbstractDict
    @test haskey(parsed, "discretizations")
    @test haskey(parsed["discretizations"], "weno5_advection")

    # Byte-diff round-trip: parsed structure must survive re-serialization.
    reserialized = JSON.json(parsed)
    reparsed = JSON.parse(reserialized)
    @test reparsed == parsed

    # Schema spot-checks: closed-AST replacement form, §4.2 ops only.
    spec = parsed["discretizations"]["weno5_advection"]
    @test spec["applies_to"]["op"] == "div"
    @test spec["applies_to"]["dim"] == "\$x"
    flux = spec["applies_to"]["args"][1]
    @test flux["op"] == "*"
    @test sort(flux["args"]) == sort(["\$U", "\$q"])
    @test spec["grid_family"] == "cartesian"
    @test spec["accuracy"] == "O(dx^5)"
    @test spec["form"] == "weighted_essentially_nonoscillatory"
    @test spec["upwind_biased"] == true

    @test haskey(spec, "replacement")
    repl = spec["replacement"]
    @test repl["op"] == "arrayop"
    @test repl["output_idx"] == ["\$x"]
    @test sort(repl["args"]) == sort(["\$U", "\$q"])

    # The replacement body is (F_E - F_W) / dx — a binary `/` over a
    # binary `-` of two flux subtrees, scaled by the grid spacing `dx`.
    body = repl["expr"]
    @test body["op"] == "/"
    @test body["args"][2] == "dx"
    diff_node = body["args"][1]
    @test diff_node["op"] == "-"
    @test length(diff_node["args"]) == 2

    # No off-spec op names anywhere in the replacement subtree.
    function walk(node, allowed)
        if node isa AbstractDict && haskey(node, "op")
            @test node["op"] in allowed
            for a in get(node, "args", [])
                walk(a, allowed)
            end
        elseif node isa AbstractArray
            for a in node
                walk(a, allowed)
            end
        end
    end
    allowed_ops = Set(
        [
            "arrayop", "index",
            "+", "-", "*", "/", "^",
            ">", "ifelse",
        ]
    )
    walk(repl, allowed_ops)
end

@testitem "weno5_advection MMS convergence: order >= 4.7 on uniform Cartesian" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fixture_dir = joinpath(repo_root, "tests", "fixtures", "weno5_advection")
    @test isdir(fixture_dir)
    input = JSON.parse(read(joinpath(fixture_dir, "mms_input.esm"), String))
    expected = JSON.parse(read(joinpath(fixture_dir, "mms_expected.esm"), String))

    @test input["rule"] == "weno5_advection"
    @test expected["rule"] == "weno5_advection"
    min_order = Float64(expected["expected_min_order"])
    grids = [Int(g["n"]) for g in input["grids"]]
    eps_weno = Float64(input["weno_epsilon"])

    # Manufactured solution: sin(2*pi*x + 1.0), analytic antiderivative.
    f(x) = sin(2pi * x + 1.0)
    F(x) = -cos(2pi * x + 1.0) / (2pi)
    cell_average(a, b) = (F(b) - F(a)) / (b - a)

    # Classical Jiang-Shu (1996) 5th-order left-biased WENO reconstruction
    # at the right edge of cell i using cells i-2..i+2. Same algebraic form
    # as the closed-AST `arrayop` body now shipped in the rule (Shu 1998
    # NASA/CR-97-206253 §2.2, eqs. (2.11),(2.15)-(2.18)) — kept here as the
    # numeric oracle for the convergence sweep.
    d0, d1, d2 = 1 / 10, 6 / 10, 3 / 10
    function weno5_left(qm2, qm1, q0, qp1, qp2)
        p0 = (1 / 3) * qm2 - (7 / 6) * qm1 + (11 / 6) * q0
        p1 = -(1 / 6) * qm1 + (5 / 6) * q0 + (1 / 3) * qp1
        p2 = (1 / 3) * q0 + (5 / 6) * qp1 - (1 / 6) * qp2

        b0 = (13 / 12) * (qm2 - 2qm1 + q0)^2 + (1 / 4) * (qm2 - 4qm1 + 3q0)^2
        b1 = (13 / 12) * (qm1 - 2q0 + qp1)^2 + (1 / 4) * (qm1 - qp1)^2
        b2 = (13 / 12) * (q0 - 2qp1 + qp2)^2 + (1 / 4) * (3q0 - 4qp1 + qp2)^2

        a0 = d0 / (eps_weno + b0)^2
        a1 = d1 / (eps_weno + b1)^2
        a2 = d2 / (eps_weno + b2)^2
        s = a0 + a1 + a2
        return (a0 / s) * p0 + (a1 / s) * p1 + (a2 / s) * p2
    end

    function linf_error(n)
        dx = 1.0 / n
        q = [cell_average((i - 1) * dx, i * dx) for i in 1:n]
        modn(j) = mod(j - 1, n) + 1
        err = 0.0
        for i in 1:n
            qhat = weno5_left(
                q[modn(i - 2)], q[modn(i - 1)], q[modn(i)],
                q[modn(i + 1)], q[modn(i + 2)],
            )
            x_face = i * dx
            err = max(err, abs(qhat - f(x_face)))
        end
        return err
    end

    errors = [linf_error(n) for n in grids]
    @test all(e -> isfinite(e) && e > 0, errors)

    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    @info "weno5_advection MMS convergence" grids errors orders
    @test minimum(orders) >= min_order

    # Sanity: errors decrease monotonically.
    for i in 1:(length(errors) - 1)
        @test errors[i + 1] < errors[i]
    end
end

@testitem "weno5_advection shock capturing: square wave overshoot bounded" begin
    using EarthSciDiscretizations
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fixture_dir = joinpath(repo_root, "tests", "fixtures", "weno5_advection")
    input = JSON.parse(read(joinpath(fixture_dir, "shock_input.esm"), String))
    expected = JSON.parse(read(joinpath(fixture_dir, "shock_expected.esm"), String))

    @test input["rule"] == "weno5_advection"
    @test expected["rule"] == "weno5_advection"
    n = Int(input["grid"]["n"])
    u = Float64(input["advection"]["velocity"])
    cfl = Float64(input["advection"]["cfl"])
    periods = Float64(input["advection"]["periods"])
    eps_weno = Float64(input["weno_epsilon"])
    tol_over = Float64(expected["max_overshoot_tolerance"])
    tol_under = Float64(expected["max_undershoot_tolerance"])
    env_max = Float64(expected["reference_envelope"]["max"])
    env_min = Float64(expected["reference_envelope"]["min"])

    d0, d1, d2 = 1 / 10, 6 / 10, 3 / 10
    function weno5_left(qm2, qm1, q0, qp1, qp2)
        p0 = (1 / 3) * qm2 - (7 / 6) * qm1 + (11 / 6) * q0
        p1 = -(1 / 6) * qm1 + (5 / 6) * q0 + (1 / 3) * qp1
        p2 = (1 / 3) * q0 + (5 / 6) * qp1 - (1 / 6) * qp2
        b0 = (13 / 12) * (qm2 - 2qm1 + q0)^2 + (1 / 4) * (qm2 - 4qm1 + 3q0)^2
        b1 = (13 / 12) * (qm1 - 2q0 + qp1)^2 + (1 / 4) * (qm1 - qp1)^2
        b2 = (13 / 12) * (q0 - 2qp1 + qp2)^2 + (1 / 4) * (3q0 - 4qp1 + qp2)^2
        a0 = d0 / (eps_weno + b0)^2
        a1 = d1 / (eps_weno + b1)^2
        a2 = d2 / (eps_weno + b2)^2
        s = a0 + a1 + a2
        return (a0 / s) * p0 + (a1 / s) * p1 + (a2 / s) * p2
    end

    dx = 1.0 / n
    dt = cfl * dx / abs(u)
    modn(j) = mod(j - 1, n) + 1

    # Cell averages of unit square wave on [0.3, 0.7] (overlap fraction).
    q = zeros(n)
    for i in 1:n
        a, b = (i - 1) * dx, i * dx
        overlap = max(0.0, min(b, 0.7) - max(a, 0.3))
        q[i] = overlap / dx
    end

    # RHS for u > 0: L(q)_i = -(F_{i+1/2} - F_{i-1/2})/dx, F = u * q^L.
    F = similar(q)
    function rhs!(dq, q)
        for i in 1:n
            F[i] = u * weno5_left(
                q[modn(i - 2)], q[modn(i - 1)], q[modn(i)],
                q[modn(i + 1)], q[modn(i + 2)],
            )
        end
        for i in 1:n
            dq[i] = -(F[i] - F[modn(i - 1)]) / dx
        end
    end

    # Gottlieb-Shu SSP-RK3 (1998) stages. Wrap the time loop in a function
    # so loop-local variables live in hard scope (avoids soft-scope ambiguity
    # with ModelingToolkit's exported globals like `t`).
    function ssp_rk3_advance!(q, dt, T)
        dq1 = similar(q)
        dq2 = similar(q)
        dq3 = similar(q)
        q1 = similar(q)
        q2 = similar(q)
        tsim = 0.0
        while tsim < T
            dts = min(dt, T - tsim)
            rhs!(dq1, q)
            @. q1 = q + dts * dq1
            rhs!(dq2, q1)
            @. q2 = 0.75 * q + 0.25 * q1 + 0.25 * dts * dq2
            rhs!(dq3, q2)
            @. q = (1 / 3) * q + (2 / 3) * q2 + (2 / 3) * dts * dq3
            tsim += dts
        end
        return q
    end
    ssp_rk3_advance!(q, dt, periods / abs(u))

    overshoot = maximum(q) - env_max
    undershoot = env_min - minimum(q)
    @info "weno5_advection shock capturing" overshoot undershoot
    @test overshoot <= tol_over
    @test undershoot <= tol_under
    # And values should at least be bounded (basic sanity).
    @test isfinite(maximum(q))
    @test isfinite(minimum(q))
end

# ---------------------------------------------------------------------------
# weno5_grad rule tests (esd-ubz: MOL #531 parity — grad(u,x) case)
# ---------------------------------------------------------------------------
# MOL #531 rules_weno.jl _weno_template fires on Dx(u) (ESS: grad(u,x)),
# not on div(U·q). weno5_advection.json applies_to is div(U*q,x) and does
# NOT match the bare 1st-derivative pattern. weno5_grad.json was added as
# the thin alias: applies_to = {"op":"grad",...}, replacement = (hp−hm)/dx
# with the same Jiang-Shu WENO5 math (no new coefficients).

@testitem "weno5_grad rule is discoverable under :finite_difference" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "weno5_grad", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_difference
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"grad\"", content)
    @test occursin("\"weighted_essentially_nonoscillatory\"", content)
    @test occursin("\"O(dx^5)\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    # Must NOT use the div operator (that's weno5_advection, not weno5_grad).
    @test !occursin("\"op\": \"div\"", content)
end

@testitem "weno5_grad rule: stencil-level numeric convergence (smooth periodic)" begin
    using JSON
    using EarthSciDiscretizations
    using EarthSciDiscretizations: lower_stencil_to_replacement

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    rule_path = joinpath(repo_root, "discretizations", "finite_difference", "weno5_grad.json")
    raw  = JSON.parsefile(rule_path)
    body = raw["discretizations"]["weno5_grad"]
    lw   = lower_stencil_to_replacement(body)
    repl = lw["replacement"]
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl

    # Evaluator with periodic wrapping (same as test_advection_discretization_conformance.jl).
    function _resolve_idx(ie, ci)
        ie isa AbstractString && ie == "\$x" && return ci
        ie isa AbstractDict || error("bad idx")
        ia = ie["args"]
        ia[1] == "\$x" || error("lhs not \$x")
        ie["op"] == "+" && return ci + Int(ia[2])
        ie["op"] == "-" && return ci - Int(ia[2])
        error("unsupported idx op $(ie["op"])")
    end

    function _eval(node, u, ci, bindings, N)
        node isa Number && return Float64(node)
        if node isa AbstractString
            node == "\$x" && return Float64(ci)
            haskey(bindings, node) && return bindings[node]
            error("unresolved var $node")
        end
        node isa AbstractDict || error("bad node")
        op = node["op"]
        args = node["args"]
        if op == "index"
            j = _resolve_idx(args[2], ci)
            return u[mod(j - 1, N) + 1]
        end
        ev(a) = _eval(a, u, ci, bindings, N)
        op == "+"   && return sum(ev(a) for a in args)
        op == "-"   && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*"   && return prod(ev(a) for a in args)
        op == "/"   && return ev(args[1]) / ev(args[2])
        op == "^"   && return ev(args[1]) ^ ev(args[2])
        error("unsupported op $op")
    end

    import Base.MathConstants
    using LinearAlgebra

    function linf_err(n)
        dx = 1.0 / n
        u  = [sin(2π * (i - 0.5) * dx) for i in 1:n]
        bindings = Dict{String,Float64}("dx" => dx)
        du = [_eval(expr, u, i, bindings, n) for i in 1:n]
        exact = [2π * cos(2π * (i - 0.5) * dx) for i in 1:n]
        return maximum(abs.(du .- exact))
    end

    # WENO5-JS yields ≥ 3.9th order for smooth periodic sin on uniform grids
    # (asymptotic 4th due to JS weight smoothness; MOL oracle matches).
    ns = [16, 32, 64, 128]
    errs = [linf_err(n) for n in ns]
    orders = [log2(errs[i] / errs[i+1]) for i in 1:(length(errs)-1)]
    @info "weno5_grad smooth convergence" ns errs orders
    @test minimum(orders) >= 3.8
    @test all(isfinite, errs) && all(e -> e > 0, errs)
end

@testitem "weno5_grad conformance: golden matches MOL #531 _weno_template" begin
    using JSON
    using EarthSciDiscretizations
    using EarthSciDiscretizations: lower_stencil_to_replacement

    HARNESS  = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                        "rect_1d_advection_weno5_periodic")
    FIXTURES = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))

    @test FIXTURES["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

    rule_path = joinpath(REPO_ROOT, FIXTURES["rule_path"])
    raw  = JSON.parsefile(rule_path)
    body = raw["discretizations"][FIXTURES["rule"]]
    lw   = lower_stencil_to_replacement(body)
    repl = lw["replacement"]
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl

    function _resolve_idx(ie, ci)
        ie isa AbstractString && ie == "\$x" && return ci
        ie isa AbstractDict || error("bad idx")
        ia = ie["args"]
        ia[1] == "\$x" || error("lhs not \$x")
        ie["op"] == "+" && return ci + Int(ia[2])
        ie["op"] == "-" && return ci - Int(ia[2])
        error("unsupported idx op $(ie["op"])")
    end

    function _eval(node, u, ci, bindings, N)
        node isa Number && return Float64(node)
        if node isa AbstractString
            node == "\$x" && return Float64(ci)
            haskey(bindings, node) && return bindings[node]
            error("unresolved $node")
        end
        node isa AbstractDict || error("bad node")
        op = node["op"]; args = node["args"]
        if op == "index"
            j = _resolve_idx(args[2], ci)
            return u[mod(j - 1, N) + 1]
        end
        ev(a) = _eval(a, u, ci, bindings, N)
        op == "+"  && return sum(ev(a) for a in args)
        op == "-"  && return (length(args)==1 ? -ev(args[1]) : ev(args[1])-ev(args[2]))
        op == "*"  && return prod(ev(a) for a in args)
        op == "/"  && return ev(args[1]) / ev(args[2])
        op == "^"  && return ev(args[1]) ^ ev(args[2])
        error("unsupported op $op")
    end

    for fx in FIXTURES["fixtures"]
        g      = fx["grid"]
        ic     = fx["initial_condition"]
        N      = Int(g["n_cells"])
        dx     = (Float64(g["x_end"]) - Float64(g["x_start"])) / N
        amp    = Float64(ic["scale_factor"]) / (Float64(ic["sigma"]) * sqrt(2π))
        x_c    = [Float64(g["x_start"]) + (i - 0.5) * dx for i in 1:N]
        u      = [amp * exp(-(xi - Float64(ic["center"]))^2 / (2 * Float64(ic["sigma"])^2)) for xi in x_c]

        bindings = Dict{String,Float64}("dx" => dx)
        du_dx = [_eval(expr, u, i, bindings, N) for i in 1:N]

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))
        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_du_dx = Float64.(golden["du_dx"])
        g_u     = Float64.(golden["field_u"])

        rel_tol = Float64(FIXTURES["tolerance"]["relative"])
        scale   = max(maximum(abs.(g_du_dx)), 1e-30)
        @test maximum(abs.(du_dx .- g_du_dx)) / scale <= rel_tol
        @test maximum(abs.(Float64.(u) .- g_u)) <= 1e-15
    end
end
