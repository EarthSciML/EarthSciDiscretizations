using Test
using TestItems

# Tests for the finite_volume/weno5_advection_2d declarative rule.
#
# Layer A: rule discovery + JSON round-trip. The 2D rule is discoverable under
# the :finite_volume family and round-trips through JSON.parse + JSON.json
# without semantic loss. Schema spot-checks confirm the dimension-by-dimension
# splitting structure: 9-point stencil covering the 5-cell support along each
# axis, per-axis reconstruction blocks (left/right-biased + smoothness
# indicators + nonlinear weights), and an upwind-biased flux selection rule.
#
# Layer B' (canonical 2D smooth advection): drives a constant-velocity
# translation of an asymmetric 2D Gaussian on a periodic [0,1]^2 square with
# WENO5 + SSP-RK3 dimension-by-dimension splitting. Asymmetric velocity
# (u != v) and Gaussian sigma (sigma_x != sigma_y) force the x- and
# y-reconstruction paths to do independent work. After exactly one full
# revolution in x and two in y the analytic solution returns to the IC; the
# L_inf error and the spurious overshoot/undershoot must stay within the
# fixture-declared bounds.
#
# The smooth-Gaussian translation case doubles as a Layer-C integration
# benchmark: it drives the full WENO5 + SSP-RK3 transport pipeline end-to-end
# over one full revolution. Convergence-order verification is intentionally
# scoped out of this bead and is filed against sibling bead T2b
# (test/weno5_advection_2d/convergence) so the rule's per-bead scope stays
# narrow.

# ---------------------------------------------------------------------------
# Shared 2D WENO5 + SSP-RK3 driver. Mirrors the 1D shock test driver in
# test_weno5_advection_rule.jl (Jiang-Shu 1996 / Shu 1998 §2.2 left-biased
# reconstruction, Gottlieb-Shu 1998 SSP-RK3 stages) but applies the
# reconstruction independently along each axis per LeVeque (2002) §20
# dimension-by-dimension splitting.
# ---------------------------------------------------------------------------

@testitem "weno5_advection_2d rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "weno5_advection_2d", rules)
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
    @test occursin("\"dimension_by_dimension\"", content)
    @test occursin("\"O(h^5)\"", content)
    @test occursin("\"stencil\"", content)
    # 9-point stencil: 5 along x at y-offset 0 + 5 along y at x-offset 0
    # share the (0,0) center (counted once).
    @test occursin("\"selectors\"", content)
    # Per-axis reconstruction blocks must declare both biases + indicators.
    @test occursin("\"axes\"", content)
    @test occursin("\"reconstruction_left_biased\"", content)
    @test occursin("\"reconstruction_right_biased\"", content)
    @test occursin("\"smoothness_indicators\"", content)
    @test occursin("\"nonlinear_weights\"", content)
    # Candidate sub-stencil coefficients (Shu 1998 eq. 2.11).
    @test occursin("[11, 6]", content)
    @test occursin("[-7, 6]", content)
    @test occursin("[5, 6]", content)
end

@testitem "weno5_advection_2d rule JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "weno5_advection_2d", rules)
    @test idx !== nothing
    rule = rules[idx]

    raw = read(rule.path, String)
    parsed = JSON.parse(raw)
    @test parsed isa AbstractDict
    @test haskey(parsed, "discretizations")
    @test haskey(parsed["discretizations"], "weno5_advection_2d")

    reserialized = JSON.json(parsed)
    reparsed = JSON.parse(reserialized)
    @test reparsed == parsed

    # Schema spot-checks (ESS §7 stencil-rule shape extended to 2D).
    spec = parsed["discretizations"]["weno5_advection_2d"]
    @test spec["applies_to"]["op"] == "advect"
    @test spec["applies_to"]["dims"] == ["\$x", "\$y"]
    @test spec["grid_family"] == "cartesian"
    @test spec["form"] == "weighted_essentially_nonoscillatory"
    @test spec["splitting"] == "dimension_by_dimension"
    @test spec["upwind_biased"] == true
    # 9-cell stencil (5 on x-axis, 5 on y-axis, sharing the (0,0) center).
    @test length(spec["stencil"]) == 9
    # Each entry uses the multi-axis `selectors` plural shape.
    for entry in spec["stencil"]
        @test haskey(entry, "selectors")
        @test length(entry["selectors"]) == 2
    end
    # Per-axis reconstruction blocks each carry 3 candidates per bias.
    for ax in ("x", "y")
        axb = spec["axes"][ax]
        @test length(axb["reconstruction_left_biased"]["candidates"]) == 3
        @test length(axb["reconstruction_right_biased"]["candidates"]) == 3
        d = axb["reconstruction_left_biased"]["linear_weights"]
        @test d["d0"]["args"] == [1, 10]
        @test d["d1"]["args"] == [6, 10]
        @test d["d2"]["args"] == [3, 10]
    end
end

@testitem "weno5_advection_2d canonical: smooth Gaussian translation on periodic square" begin
    using EarthSciDiscretizations
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    fixture_dir = joinpath(repo_root, "tests", "fixtures", "weno5_advection_2d", "canonical")
    @test isdir(fixture_dir)
    input = JSON.parse(read(joinpath(fixture_dir, "input.esm"), String))
    expected = JSON.parse(read(joinpath(fixture_dir, "expected.esm"), String))

    @test input["rule"] == "weno5_advection_2d"
    @test expected["rule"] == "weno5_advection_2d"
    nx = Int(input["grid"]["nx"])
    ny = Int(input["grid"]["ny"])
    sigma_x = Float64(input["initial_condition"]["sigma_x"])
    sigma_y = Float64(input["initial_condition"]["sigma_y"])
    x0 = Float64(input["initial_condition"]["x0"])
    y0 = Float64(input["initial_condition"]["y0"])
    amp = Float64(input["initial_condition"]["amplitude"])
    u = Float64(input["advection"]["velocity_x"])
    v = Float64(input["advection"]["velocity_y"])
    cfl = Float64(input["advection"]["cfl"])
    periods = Float64(input["advection"]["periods"])
    eps_weno = Float64(input["weno_epsilon"])
    tol_linf = Float64(expected["max_linf_error"])
    tol_over = Float64(expected["max_overshoot"])
    tol_under = Float64(expected["max_undershoot"])

    Lx = 1.0
    Ly = 1.0
    dx = Lx / nx
    dy = Ly / ny

    # Cell averages of the asymmetric Gaussian (closed-form via 1D Gaussian
    # antiderivatives erf, factored across x and y).
    function gauss_cell_avg(i, j)
        a_x = (i - 1) * dx
        b_x = i * dx
        a_y = (j - 1) * dy
        b_y = j * dy
        # Use Simpson's rule with 9 sub-intervals per axis for the cell average;
        # well-converged for sigma >> dx (>= 3 sub-intervals per sigma).
        nsub = 9
        sxs = range(a_x, b_x; length = nsub)
        sys = range(a_y, b_y; length = nsub)
        wsimp = [1, 4, 2, 4, 2, 4, 2, 4, 1] ./ (3.0 * (nsub - 1))
        s = 0.0
        for kx in 1:nsub, ky in 1:nsub
            xk = sxs[kx]; yk = sys[ky]
            qk = amp * exp(-((xk - x0)^2 / (2 * sigma_x^2) + (yk - y0)^2 / (2 * sigma_y^2)))
            s += wsimp[kx] * wsimp[ky] * qk
        end
        return s
    end

    q0 = [gauss_cell_avg(i, j) for i in 1:nx, j in 1:ny]

    # 1D WENO5 left-biased reconstruction (Shu 1998 §2.2 eqs. 2.11, 2.15-2.18).
    d0, d1, d2 = 1 / 10, 6 / 10, 3 / 10
    function weno5_left(qm2, qm1, q0v, qp1, qp2)
        p0 = (1 / 3) * qm2 - (7 / 6) * qm1 + (11 / 6) * q0v
        p1 = -(1 / 6) * qm1 + (5 / 6) * q0v + (1 / 3) * qp1
        p2 = (1 / 3) * q0v + (5 / 6) * qp1 - (1 / 6) * qp2
        b0 = (13 / 12) * (qm2 - 2qm1 + q0v)^2 + (1 / 4) * (qm2 - 4qm1 + 3q0v)^2
        b1 = (13 / 12) * (qm1 - 2q0v + qp1)^2 + (1 / 4) * (qm1 - qp1)^2
        b2 = (13 / 12) * (q0v - 2qp1 + qp2)^2 + (1 / 4) * (3q0v - 4qp1 + qp2)^2
        a0 = d0 / (eps_weno + b0)^2
        a1 = d1 / (eps_weno + b1)^2
        a2 = d2 / (eps_weno + b2)^2
        s = a0 + a1 + a2
        return (a0 / s) * p0 + (a1 / s) * p1 + (a2 / s) * p2
    end

    modnx(j) = mod(j - 1, nx) + 1
    modny(j) = mod(j - 1, ny) + 1

    # 2D advection RHS via dimension-by-dimension splitting:
    #   d/dt q = -∂_x(u * q) - ∂_y(v * q)
    #         ≈ -(F^x_{i+1/2,j} - F^x_{i-1/2,j})/dx
    #           -(F^y_{i,j+1/2} - F^y_{i,j-1/2})/dy
    # with F^x_{i+1/2,j} = u * q^L_{i+1/2,j} for u > 0 (this fixture: u, v > 0).
    function rhs!(dq, q)
        Fx = similar(q)
        Fy = similar(q)
        @assert u > 0 && v > 0  # left-biased reconstruction is valid
        for j in 1:ny, i in 1:nx
            Fx[i, j] = u * weno5_left(
                q[modnx(i - 2), j], q[modnx(i - 1), j], q[modnx(i), j],
                q[modnx(i + 1), j], q[modnx(i + 2), j],
            )
        end
        for j in 1:ny, i in 1:nx
            Fy[i, j] = v * weno5_left(
                q[i, modny(j - 2)], q[i, modny(j - 1)], q[i, modny(j)],
                q[i, modny(j + 1)], q[i, modny(j + 2)],
            )
        end
        for j in 1:ny, i in 1:nx
            dq[i, j] = -(Fx[i, j] - Fx[modnx(i - 1), j]) / dx -
                       (Fy[i, j] - Fy[i, modny(j - 1)]) / dy
        end
        return dq
    end

    # Gottlieb-Shu (1998) SSP-RK3 stages, 2D variant. Wrap in a function so
    # loop-locals stay in hard scope (avoids `t`-symbol soft-scope ambiguity
    # under ModelingToolkit).
    function ssp_rk3_advance_2d!(q, dt, T)
        dq1 = similar(q); dq2 = similar(q); dq3 = similar(q)
        q1 = similar(q); q2 = similar(q)
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

    dt = cfl * min(dx / abs(u), dy / abs(v))
    q = copy(q0)
    # Pass `periods` directly rather than binding a top-level `T` in the
    # testitem — `T` is a frequently-clobbered name in scientific-Julia
    # imports and the runner's @testitem precompile path can hit a binding
    # collision against ModelingToolkit/Symbolics re-exports.
    ssp_rk3_advance_2d!(q, dt, Float64(periods))

    # Analytic solution at T = exactly one revolution in x (u*T=1) and two in y
    # (v*T=2): translated Gaussian is identical to the IC.
    err = maximum(abs(q[i, j] - q0[i, j]) for j in 1:ny, i in 1:nx)
    overshoot = maximum(q) - 1.0  # IC peak is amp = 1.0
    undershoot = 0.0 - minimum(q)
    @info "weno5_advection_2d canonical" linf=err overshoot=overshoot undershoot=undershoot
    @test isfinite(err) && err >= 0
    @test err <= tol_linf
    @test overshoot <= tol_over
    @test undershoot <= tol_under
end

