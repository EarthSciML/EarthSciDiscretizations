@testsnippet AdvectionDiscretizationConformanceSetup begin
    using Test
    using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
    using JSON
end

# ---------------------------------------------------------------------------
# Helpers: shared canonical-pipeline evaluator for 1D periodic advection
# conformance. Duplicates the logic in
# tests/conformance/discretization/*/regenerate_golden.jl intentionally — the
# conformance harness must be independent of the generator so that regressions
# in ESD are visible (if the generator and test share code, a bug in the shared
# code can mask a regression).
# ---------------------------------------------------------------------------

@testsnippet AdvectionConformanceHelpers begin
    using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
    using JSON

    function _has_index(node)::Bool
        node isa AbstractDict || return false
        node["op"] == "index" && return true
        args = get(node, "args", Any[])
        return any(_has_index(a) for a in args)
    end

    function _eval_replacement(
            node,
            u::Vector{Float64},
            cell_idx::Int,
            bindings::Dict{String, Float64},
            N::Int,
    )::Float64
        node isa Number && return Float64(node)
        if node isa AbstractString
            s = String(node)
            s == "\$x" && return Float64(cell_idx)
            haskey(bindings, s) && return bindings[s]
            error("_eval_replacement: unresolved variable '$s'")
        end
        node isa AbstractDict || error("_eval_replacement: unexpected node type")
        op   = String(node["op"])
        args = node["args"]
        if op == "index"
            function _resolve_idx(ie)
                ie isa AbstractString && String(ie) == "\$x" && return cell_idx
                ie isa AbstractDict || error("unsupported index: $ie")
                ia = ie["args"]
                String(ia[1]) == "\$x" || error("index lhs not '\$x': $(ia[1])")
                op2 = String(ie["op"])
                op2 == "+" && return cell_idx + Int(ia[2])
                op2 == "-" && return cell_idx - Int(ia[2])
                error("unsupported index op '$op2'")
            end
            j = _resolve_idx(args[2])
            return u[mod(j - 1, N) + 1]
        end
        !_has_index(node) && return eval_coeff(node, bindings)
        ev(a) = _eval_replacement(a, u, cell_idx, bindings, N)
        op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
        op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _apply_rule_1d_periodic(rule_path::String, rule_name::String,
                                     u::Vector{Float64}, dx::Float64)
        raw  = JSON.parsefile(rule_path)
        body = raw["discretizations"][rule_name]
        lw   = lower_stencil_to_replacement(body)
        repl = lw["replacement"]
        expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
               repl["expr"] : repl
        N    = length(u)
        bindings = Dict{String, Float64}("dx" => dx)
        return Float64[_eval_replacement(expr, u, i, bindings, N) for i in 1:N]
    end

    function _gaussian_field(n::Int, x_start::Float64, x_end::Float64,
                             center::Float64, sigma::Float64, scale::Float64)
        dx  = (x_end - x_start) / n
        amp = scale / (sigma * sqrt(2π))
        return Float64[amp * exp(-(x_start + (i - 0.5) * dx - center)^2 /
                                 (2 * sigma^2)) for i in 1:n]
    end
end

@testsnippet NonuniformAdvectionHelpers begin
    using EarthSciDiscretizations: lower_stencil_to_replacement
    using JSON

    # Per-cell evaluator for rules with index("dx", "$x") coefficients.
    function _resolve_idx_nu(ie, cell_idx::Int)::Int
        ie isa AbstractString && String(ie) == "\$x" && return cell_idx
        ie isa AbstractDict || error("unsupported index: $ie")
        ia = ie["args"]
        String(ia[1]) == "\$x" || error("index LHS must be '\$x'")
        op2 = String(ie["op"])
        op2 == "+" && return cell_idx + Int(ia[2])
        op2 == "-" && return cell_idx - Int(ia[2])
        error("unsupported index op '$op2'")
    end

    function _eval_nonuniform(node, u::Vector{Float64}, dx_arr::Vector{Float64},
                              cell_idx::Int, N::Int)::Float64
        node isa Number && return Float64(node)
        if node isa AbstractString
            s = String(node)
            s == "\$x" && return Float64(cell_idx)
            error("_eval_nonuniform: unresolved variable '$s'")
        end
        node isa AbstractDict || error("unexpected node type: $(typeof(node))")
        op   = String(node["op"])
        args = node["args"]
        if op == "index"
            arr_name = String(args[1])
            j = _resolve_idx_nu(args[2], cell_idx)
            arr_name == "\$u" && return u[mod(j - 1, N) + 1]
            arr_name == "dx"  && return dx_arr[mod(j - 1, N) + 1]
            error("unknown array '$arr_name'")
        end
        ev(a) = _eval_nonuniform(a, u, dx_arr, cell_idx, N)
        op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
        op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*" && return prod(ev(a) for a in args)
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _apply_rule_1d_nonuniform_periodic(rule_path::String, rule_name::String,
                                                u::Vector{Float64}, dx_arr::Vector{Float64})
        raw  = JSON.parsefile(rule_path)
        body = raw["discretizations"][rule_name]
        lw   = lower_stencil_to_replacement(body)
        repl = lw["replacement"]
        expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
               repl["expr"] : repl
        N = length(u)
        return Float64[_eval_nonuniform(expr, u, dx_arr, i, N) for i in 1:N]
    end

    function _stretched_dx(n::Int, x_start::Float64, x_end::Float64,
                           amplitude::Float64)::Vector{Float64}
        L = x_end - x_start
        h = L / n
        dx = Float64[h * (1.0 + amplitude * sin(2π * (i - 1) / n)) for i in 1:n]
        @assert abs(sum(dx) - L) < 1e-10
        return dx
    end

    function _gaussian_field_nonuniform(dx_arr::Vector{Float64}, x_start::Float64,
                                        center::Float64, sigma::Float64,
                                        scale::Float64)::Vector{Float64}
        N = length(dx_arr)
        x_c = Vector{Float64}(undef, N)
        pos = x_start
        for i in 1:N
            x_c[i] = pos + dx_arr[i] / 2
            pos += dx_arr[i]
        end
        amp = scale / (sigma * sqrt(2π))
        return Float64[amp * exp(-(x_c[i] - center)^2 / (2 * sigma^2)) for i in 1:N]
    end
end

# ---------------------------------------------------------------------------
# Centered 2nd-order advection conformance
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_1d_advection_centered_periodic" setup=[AdvectionConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                       "rect_1d_advection_centered_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    # Verify the SHA header is present in fixtures so capture provenance is self-documenting.
    @test get(FIXTURES, "_mol531_sha", nothing) == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

    rule_path = joinpath(REPO_ROOT, FIXTURES["rule_path"])

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        N   = Int(g["n_cells"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / N
        u   = _gaussian_field(N, Float64(g["x_start"]), Float64(g["x_end"]),
                              Float64(ic["center"]), Float64(ic["sigma"]),
                              Float64(ic["scale_factor"]))

        du_dx = _apply_rule_1d_periodic(rule_path, FIXTURES["rule"], u, dx)

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        # Verify SHA header in golden.
        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_du_dx = Float64.(golden["du_dx"])
        g_u     = Float64.(golden["field_u"])
        @test length(du_dx) == length(g_du_dx)
        @test length(u)     == length(g_u)

        @testset "fixture=$(fx["name"])" begin
            for i in 1:N
                scale = max(1.0, abs(du_dx[i]), abs(g_du_dx[i]))
                @test abs(du_dx[i] - g_du_dx[i]) <= REL_TOL * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 1st-order upwind advection conformance
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_1d_advection_upwind_periodic" setup=[AdvectionConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                       "rect_1d_advection_upwind_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    rule_path = joinpath(REPO_ROOT, FIXTURES["rule_path"])

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        N   = Int(g["n_cells"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / N
        u   = _gaussian_field(N, Float64(g["x_start"]), Float64(g["x_end"]),
                              Float64(ic["center"]), Float64(ic["sigma"]),
                              Float64(ic["scale_factor"]))

        du_dx = _apply_rule_1d_periodic(rule_path, FIXTURES["rule"], u, dx)

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_du_dx = Float64.(golden["du_dx"])
        @test length(du_dx) == length(g_du_dx)

        @testset "fixture=$(fx["name"])" begin
            for i in 1:N
                scale = max(1.0, abs(du_dx[i]), abs(g_du_dx[i]))
                @test abs(du_dx[i] - g_du_dx[i]) <= REL_TOL * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 1st-order upwind non-uniform advection conformance (esd-02z)
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_1d_advection_upwind_nonuniform_periodic" setup=[NonuniformAdvectionHelpers] tags=[:conformance, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                       "rect_1d_advection_upwind_nonuniform_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    rule_path = joinpath(REPO_ROOT, FIXTURES["rule_path"])

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        N   = Int(g["n_cells"])
        x_start = Float64(g["x_start"])
        x_end   = Float64(g["x_end"])
        amp = Float64(g["stretching_amplitude"])

        dx_arr = _stretched_dx(N, x_start, x_end, amp)
        u      = _gaussian_field_nonuniform(dx_arr, x_start,
                                            Float64(ic["center"]),
                                            Float64(ic["sigma"]),
                                            Float64(ic["scale_factor"]))

        du_dx = _apply_rule_1d_nonuniform_periodic(rule_path, FIXTURES["rule"], u, dx_arr)

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test get(golden, "_captured_by", nothing) == "esd-02z"

        g_du_dx = Float64.(golden["du_dx"])
        @test length(du_dx) == length(g_du_dx)

        @testset "fixture=$(fx["name"])" begin
            for i in 1:N
                scale = max(1.0, abs(du_dx[i]), abs(g_du_dx[i]))
                @test abs(du_dx[i] - g_du_dx[i]) <= REL_TOL * scale
            end
        end
    end
end
