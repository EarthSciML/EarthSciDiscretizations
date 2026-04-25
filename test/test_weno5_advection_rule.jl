using Test
using TestItems

# Tests for the finite_volume/weno5_advection declarative rule.
#
# Layer A: rule discovery + JSON round-trip. The rule is discoverable under
# the :finite_volume family and its JSON file round-trips through
# JSON.parse + JSON.json without semantic loss.
#
# Layer B: MMS convergence. The classical Jiang-Shu (1996) 5th-order WENO
# reconstruction on a uniform 1D Cartesian grid is theoretically O(dx^5)
# in L_inf on smooth data away from critical points. Measured on the
# left-biased right-edge reconstruction q_{i+1/2}^L of a phase-shifted
# sine at n = 32, 64, 128, 256; minimum observed order must be ≥ 4.7.
#
# Layer B (shock): Linear advection of a unit square wave by u = 1 for one
# full period on a periodic domain with WENO5 + SSP-RK3. The exact solution
# is the initial condition; max overshoot and undershoot at cell faces
# must stay bounded by 0.05 (qualitative shock-capturing check).

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
    @test occursin("\"op\": \"advect\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"weighted_essentially_nonoscillatory\"", content)
    @test occursin("\"O(dx^5)\"", content)
    @test occursin("\"stencil\"", content)
    # 5-point support stencil i-2..i+2.
    for off in (-2, -1, 0, 1, 2)
        @test occursin("\"offset\": $off", content)
    end
    # Candidate sub-stencil coefficients (Shu 1998 eq. 2.11) must appear.
    @test occursin("\"reconstruction_left_biased\"", content)
    @test occursin("[11, 6]", content)
    @test occursin("[-7, 6]", content)
    @test occursin("[5, 6]", content)
    @test occursin("\"smoothness_indicators\"", content)
    @test occursin("\"nonlinear_weights\"", content)
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

    # Schema spot-checks (ESS §7 stencil-rule shape).
    spec = parsed["discretizations"]["weno5_advection"]
    @test spec["applies_to"]["op"] == "advect"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["grid_family"] == "cartesian"
    @test spec["accuracy"] == "O(dx^5)"
    @test spec["form"] == "weighted_essentially_nonoscillatory"
    @test spec["upwind_biased"] == true
    @test length(spec["stencil"]) == 5
    offsets = sort([s["selector"]["offset"] for s in spec["stencil"]])
    @test offsets == [-2, -1, 0, 1, 2]
    @test length(spec["reconstruction_left_biased"]["candidates"]) == 3
    @test length(spec["reconstruction_right_biased"]["candidates"]) == 3
    # Linear weights (Shu 1998 eq. 2.15): d0=1/10, d1=6/10, d2=3/10, sum=1.
    d = spec["reconstruction_left_biased"]["linear_weights"]
    @test d["d0"]["args"] == [1, 10]
    @test d["d1"]["args"] == [6, 10]
    @test d["d2"]["args"] == [3, 10]
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
    # at the right edge of cell i using cells i-2..i+2. Weights per
    # Shu (1998) NASA/CR-97-206253 §2.2, eqs. (2.11),(2.15)-(2.18).
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
