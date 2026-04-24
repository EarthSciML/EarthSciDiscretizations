using Test
using TestItems

# Tests for the finite_volume/ppm_reconstruction declarative rule.
#
# Layer A (rule discovery + byte-diff round-trip): the rule is discoverable
# under the :finite_volume family and its JSON file round-trips through
# JSON.parse + JSON.json without semantic loss.
#
# Layer B (MMS convergence): the unlimited PPM reconstruction (CW84 eqs.
# 1.6-1.10) on a uniform 1D Cartesian grid is theoretically O(dx^3) in
# L_inf. We exercise it on a smooth periodic manufactured solution at
# n = 16, 32, 64, 128 cells and assert minimum measured order >= 2.7.

@testitem "ppm_reconstruction rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"reconstruct\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"piecewise_parabolic\"", content)
    @test occursin("\"O(dx^3)\"", content)
    @test occursin("\"stencil\"", content)
    # 5-point support stencil i-2..i+2.
    for off in (-2, -1, 0, 1, 2)
        @test occursin("\"offset\": $off", content)
    end
    # 4th-order edge interpolation coefficients (CW84 eq. 1.6).
    @test occursin("\"edge_value_stencil\"", content)
    @test occursin("[-1, 12]", content)
    @test occursin("[7, 12]", content)
end

@testitem "ppm_reconstruction rule JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction", rules)
    @test idx !== nothing
    rule = rules[idx]

    raw = read(rule.path, String)
    parsed = JSON.parse(raw)
    @test parsed isa AbstractDict
    @test haskey(parsed, "discretizations")
    @test haskey(parsed["discretizations"], "ppm_reconstruction")

    # Byte-diff round-trip: re-parse the serialized form and check the
    # parsed structure matches the original exactly. (Whitespace differs
    # between JSON.json and the hand-formatted source file, so we compare
    # parsed structures rather than raw strings.)
    reserialized = JSON.json(parsed)
    reparsed = JSON.parse(reserialized)
    @test reparsed == parsed

    # Schema spot-checks (ESS §7 stencil-rule shape).
    spec = parsed["discretizations"]["ppm_reconstruction"]
    @test spec["applies_to"]["op"] == "reconstruct"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["grid_family"] == "cartesian"
    @test spec["accuracy"] == "O(dx^3)"
    @test spec["form"] == "piecewise_parabolic"
    @test spec["limiter"] == "none"
    @test length(spec["stencil"]) == 5
    offsets = sort([s["selector"]["offset"] for s in spec["stencil"]])
    @test offsets == [-2, -1, 0, 1, 2]
    edge = spec["edge_value_stencil"]
    @test length(edge["stencil"]) == 4
end

@testitem "ppm_reconstruction MMS convergence: order >= 2.7 on uniform Cartesian" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    # --- Test setup ---------------------------------------------------------
    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fixture_dir = joinpath(repo_root, "tests", "fixtures", "ppm_reconstruction")
    @test isdir(fixture_dir)
    input = JSON.parse(read(joinpath(fixture_dir, "input.esm"), String))
    expected = JSON.parse(read(joinpath(fixture_dir, "expected.esm"), String))

    @test input["rule"] == "ppm_reconstruction"
    @test expected["rule"] == "ppm_reconstruction"
    min_order = Float64(expected["expected_min_order"])
    samples_per_cell = Int(input["samples_per_cell"])
    grids = [Int(g["n"]) for g in input["grids"]]

    # Manufactured solution and its analytical antiderivative for exact
    # cell averages: f(x) = sin(2*pi*x) + 0.3*cos(4*pi*x).
    f(x) = sin(2pi * x) + 0.3 * cos(4pi * x)
    F(x) = -cos(2pi * x) / (2pi) + 0.3 * sin(4pi * x) / (4pi)
    cell_average(a, b) = (F(b) - F(a)) / (b - a)

    # 4th-order edge interpolation (CW84 eq. 1.6): q_{i+1/2} from q_{i-1..i+2}.
    edge_value(qm1, q0, qp1, qp2) = (-qm1 + 7 * q0 + 7 * qp1 - qp2) / 12

    # Parabolic reconstruction (CW84 eqs. 1.5,1.7,1.10): given a_L, a_R, q_i,
    # evaluate a(xi) for xi in [0,1].
    function parabola(a_L, a_R, q_i, xi)
        da = a_R - a_L
        a6 = 6 * (q_i - 0.5 * (a_L + a_R))
        return a_L + xi * (da + a6 * (1 - xi))
    end

    # --- Per-grid L_inf error ----------------------------------------------
    function linf_error(n)
        dx = 1.0 / n
        # Cell averages from analytical integration (no quadrature error).
        q = [cell_average((i - 1) * dx, i * dx) for i in 1:n]
        # Right-edge value q_{i+1/2}; index i=0..n-1 maps to cell i+1 in 1-based.
        modn(j) = mod(j - 1, n) + 1   # periodic 1-based index helper
        q_edge = zeros(n)
        for i in 1:n
            q_edge[i] = edge_value(q[modn(i - 1)], q[modn(i)], q[modn(i + 1)], q[modn(i + 2)])
        end
        # Sample parabola at interior sub-points within each cell.
        err = 0.0
        for i in 1:n
            a_L = q_edge[modn(i - 1)]   # left edge of cell i = right edge of cell i-1
            a_R = q_edge[i]             # right edge of cell i
            q_i = q[i]
            for k in 1:samples_per_cell
                xi = (k - 0.5) / samples_per_cell
                x = (i - 1) * dx + xi * dx
                err = max(err, abs(parabola(a_L, a_R, q_i, xi) - f(x)))
            end
        end
        return err
    end

    errors = [linf_error(n) for n in grids]
    @test all(e -> isfinite(e) && e > 0, errors)

    # --- Convergence orders -------------------------------------------------
    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    @info "ppm_reconstruction MMS convergence" grids errors orders
    @test minimum(orders) >= min_order

    # Sanity: errors decrease monotonically.
    for i in 1:(length(errors) - 1)
        @test errors[i + 1] < errors[i]
    end
end
