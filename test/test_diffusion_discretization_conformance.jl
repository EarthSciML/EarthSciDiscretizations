@testsnippet DiffusionDiscretizationConformanceSetup begin
    using Test
    using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
    using JSON
end

# ---------------------------------------------------------------------------
# Helpers: 2D periodic evaluator for arakawa-lowered stencils.
# Independent of the generator in regenerate_golden.jl — regressions in ESD
# are visible only when test and generator don't share code.
# ---------------------------------------------------------------------------

@testsnippet Diffusion2DConformanceHelpers begin
    using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
    using JSON

    function _has_index_2d(node)::Bool
        node isa AbstractDict || return false
        node["op"] == "index" && return true
        args = get(node, "args", Any[])
        return any(_has_index_2d(a) for a in args)
    end

    function _resolve_axis_idx_2d(ie, xi::Int, yi::Int)::Int
        if ie isa AbstractString
            s = String(ie)
            s == "\$x" && return xi
            s == "\$y" && return yi
            error("unknown axis variable '$s'")
        end
        ie isa AbstractDict || error("unexpected type $(typeof(ie))")
        ia  = ie["args"]
        op2 = String(ie["op"])
        lhs = String(ia[1])
        off = Int(ia[2])
        base = lhs == "\$x" ? xi : (lhs == "\$y" ? yi : error("unknown lhs '$lhs'"))
        op2 == "+" && return base + off
        op2 == "-" && return base - off
        error("unsupported offset op '$op2'")
    end

    function _eval_replacement_2d(
            node,
            u::Matrix{Float64},
            xi::Int, yi::Int,
            bindings::Dict{String, Float64},
            Nx::Int, Ny::Int,
    )::Float64
        node isa Number && return Float64(node)
        if node isa AbstractString
            s = String(node)
            s == "\$x" && return Float64(xi)
            s == "\$y" && return Float64(yi)
            haskey(bindings, s) && return bindings[s]
            error("unresolved variable '$s'")
        end
        node isa AbstractDict || error("unexpected node type $(typeof(node))")
        op   = String(node["op"])
        args = node["args"]
        if op == "index"
            length(args) == 3 || error("index expects 3 args, got $(length(args))")
            jx = _resolve_axis_idx_2d(args[2], xi, yi)
            jy = _resolve_axis_idx_2d(args[3], xi, yi)
            return u[mod(jx - 1, Nx) + 1, mod(jy - 1, Ny) + 1]
        end
        !_has_index_2d(node) && return eval_coeff(node, bindings)
        ev(a) = _eval_replacement_2d(a, u, xi, yi, bindings, Nx, Ny)
        op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
        op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _gaussian_2d_field(nx::Int, ny::Int,
                                x_start::Float64, x_end::Float64,
                                y_start::Float64, y_end::Float64,
                                cx::Float64, cy::Float64,
                                sigma::Float64, scale::Float64)
        dx  = (x_end - x_start) / nx
        dy  = (y_end - y_start) / ny
        amp = scale / (2π * sigma^2)
        u   = Matrix{Float64}(undef, nx, ny)
        for ix in 1:nx
            x = x_start + (ix - 0.5) * dx
            for iy in 1:ny
                y = y_start + (iy - 0.5) * dy
                u[ix, iy] = amp * exp(-((x - cx)^2 + (y - cy)^2) / (2 * sigma^2))
            end
        end
        return u
    end
end

# ---------------------------------------------------------------------------
# Helper: 1D periodic evaluator for single-field 1D rules (shared with
# Burgers diffusion case below).
# ---------------------------------------------------------------------------

@testsnippet Diffusion1DConformanceHelpers begin
    using EarthSciDiscretizations: lower_stencil_to_replacement, eval_coeff
    using JSON

    function _has_index_1d(node)::Bool
        node isa AbstractDict || return false
        node["op"] == "index" && return true
        args = get(node, "args", Any[])
        return any(_has_index_1d(a) for a in args)
    end

    function _eval_replacement_1d(
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
            error("unresolved variable '$s'")
        end
        node isa AbstractDict || error("unexpected node type $(typeof(node))")
        op   = String(node["op"])
        args = node["args"]
        if op == "index"
            function _resolve_idx(ie)
                ie isa AbstractString && String(ie) == "\$x" && return cell_idx
                ie isa AbstractDict || error("unsupported index: $ie")
                ia  = ie["args"]
                String(ia[1]) == "\$x" || error("index lhs not '\$x': $(ia[1])")
                op2 = String(ie["op"])
                op2 == "+" && return cell_idx + Int(ia[2])
                op2 == "-" && return cell_idx - Int(ia[2])
                error("unsupported index op '$op2'")
            end
            j = _resolve_idx(args[2])
            return u[mod(j - 1, N) + 1]
        end
        !_has_index_1d(node) && return eval_coeff(node, bindings)
        ev(a) = _eval_replacement_1d(a, u, cell_idx, bindings, N)
        op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
        op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _gaussian_field_1d(n::Int, x_start::Float64, x_end::Float64,
                                center::Float64, sigma::Float64, scale::Float64)
        dx  = (x_end - x_start) / n
        amp = scale / (sigma * sqrt(2π))
        return Float64[amp * exp(-(x_start + (i - 0.5) * dx - center)^2 / (2 * sigma^2))
                       for i in 1:n]
    end
end

# ---------------------------------------------------------------------------
# Helper: two-field 1D evaluator for nonlinear_laplacian_uniform.
# ---------------------------------------------------------------------------

@testsnippet NonlinearDiffusion1DConformanceHelpers begin
    using EarthSciDiscretizations: eval_coeff
    using JSON

    function _has_index_nl(node)::Bool
        node isa AbstractDict || return false
        node["op"] == "index" && return true
        args = get(node, "args", Any[])
        return any(_has_index_nl(a) for a in args)
    end

    function _eval_replacement_nl(
            node,
            u::Vector{Float64},
            f::Vector{Float64},
            cell_idx::Int,
            bindings::Dict{String, Float64},
            N::Int,
    )::Float64
        node isa Number && return Float64(node)
        if node isa AbstractString
            s = String(node)
            s == "\$x" && return Float64(cell_idx)
            haskey(bindings, s) && return bindings[s]
            error("unresolved variable '$s'")
        end
        node isa AbstractDict || error("unexpected node type $(typeof(node))")
        op   = String(node["op"])
        args = node["args"]
        if op == "index"
            field_name = String(args[1])
            field = field_name == "\$u" ? u :
                    field_name == "\$f" ? f :
                    error("unknown field '$field_name'")
            function _resolve_idx(ie)
                ie isa AbstractString && String(ie) == "\$x" && return cell_idx
                ie isa AbstractDict || error("unsupported index: $ie")
                ia  = ie["args"]
                String(ia[1]) == "\$x" || error("index lhs not '\$x': $(ia[1])")
                op2 = String(ie["op"])
                op2 == "+" && return cell_idx + Int(ia[2])
                op2 == "-" && return cell_idx - Int(ia[2])
                error("unsupported index op '$op2'")
            end
            j = _resolve_idx(args[2])
            return field[mod(j - 1, N) + 1]
        end
        !_has_index_nl(node) && return eval_coeff(node, bindings)
        ev(a) = _eval_replacement_nl(a, u, f, cell_idx, bindings, N)
        op == "+" && return (length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args))
        op == "-" && return (length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2]))
        op == "*" && return (length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args))
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _sinusoidal_plus_const_field(n::Int, x_start::Float64, x_end::Float64,
                                          amplitude::Float64, constant::Float64)
        dx = (x_end - x_start) / n
        return Float64[constant + amplitude * sin(π * (x_start + (i - 0.5) * dx))
                       for i in 1:n]
    end
end

# ---------------------------------------------------------------------------
# 2D diffusion (5-point Laplacian) conformance
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_2d_diffusion_5point_periodic" setup=[Diffusion2DConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS   = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                         "rect_2d_diffusion_5point_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    @test get(FIXTURES, "_mol531_sha", nothing) == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

    raw_rule  = JSON.parsefile(joinpath(REPO_ROOT, FIXTURES["rule_path"]))
    rule_name = FIXTURES["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw  = rule_lowered["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        Nx  = Int(g["n_cells_x"])
        Ny  = Int(g["n_cells_y"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / Nx
        dy  = (Float64(g["y_end"]) - Float64(g["y_start"])) / Ny

        u = _gaussian_2d_field(
            Nx, Ny,
            Float64(g["x_start"]), Float64(g["x_end"]),
            Float64(g["y_start"]), Float64(g["y_end"]),
            Float64(ic["center_x"]), Float64(ic["center_y"]),
            Float64(ic["sigma"]), Float64(ic["scale_factor"]),
        )

        bindings = Dict{String, Float64}("dx" => dx, "dy" => dy)

        lap_u = Matrix{Float64}(undef, Nx, Ny)
        for iy in 1:Ny, ix in 1:Nx
            lap_u[ix, iy] = _eval_replacement_2d(replacement, u, ix, iy, bindings, Nx, Ny)
        end

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_lap = golden["laplacian_u"]
        @test length(g_lap) == Nx
        @test length(g_lap[1]) == Ny

        @testset "fixture=$(fx["name"])" begin
            for ix in 1:Nx, iy in 1:Ny
                v     = lap_u[ix, iy]
                g_v   = Float64(g_lap[ix][iy])
                scale = max(1.0, abs(v), abs(g_v))
                @test abs(v - g_v) <= REL_TOL * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 1D Burgers diffusion (d2u/dx2, centered 2nd-order) conformance
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_1d_burgers_diffusion_periodic" setup=[Diffusion1DConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS   = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                         "rect_1d_burgers_diffusion_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    @test get(FIXTURES, "_mol531_sha", nothing) == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

    raw_rule  = JSON.parsefile(joinpath(REPO_ROOT, FIXTURES["rule_path"]))
    rule_name = FIXTURES["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw  = rule_lowered["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        ic  = fx["initial_condition"]
        N   = Int(g["n_cells"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / N

        u = _gaussian_field_1d(N, Float64(g["x_start"]), Float64(g["x_end"]),
                               Float64(ic["center"]), Float64(ic["sigma"]),
                               Float64(ic["scale_factor"]))

        bindings  = Dict{String, Float64}("dx" => dx)
        d2u_dx2   = [_eval_replacement_1d(replacement, u, i, bindings, N) for i in 1:N]

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_d2u = Float64.(golden["d2u_dx2"])
        @test length(d2u_dx2) == length(g_d2u)

        @testset "fixture=$(fx["name"])" begin
            for i in 1:N
                scale = max(1.0, abs(d2u_dx2[i]), abs(g_d2u[i]))
                @test abs(d2u_dx2[i] - g_d2u[i]) <= REL_TOL * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 1D nonlinear diffusion (nonlinear_laplacian_uniform) conformance
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_1d_nonlinear_diffusion_periodic" setup=[NonlinearDiffusion1DConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS   = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                         "rect_1d_nonlinear_diffusion_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    @test get(FIXTURES, "_mol531_sha", nothing) == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

    raw_rule  = JSON.parsefile(joinpath(REPO_ROOT, FIXTURES["rule_path"]))
    rule_name = FIXTURES["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    repl_raw  = rule_body["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    for fx in FIXTURES["fixtures"]
        g    = fx["grid"]
        ic   = fx["initial_condition"]
        N    = Int(g["n_cells"])
        dx   = (Float64(g["x_end"]) - Float64(g["x_start"])) / N

        u_ic = ic["u"]
        u = _sinusoidal_plus_const_field(
            N, Float64(g["x_start"]), Float64(g["x_end"]),
            Float64(u_ic["amplitude"]), Float64(u_ic["constant"]),
        )
        f = 1.0 ./ u

        bindings  = Dict{String, Float64}("dx" => dx)
        nl_lap_u  = [_eval_replacement_nl(replacement, u, f, i, bindings, N) for i in 1:N]

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test golden["_mol531_sha"] == "35cc9143dc553ac7d3619738bd77b250c1ed162f"

        g_nl = Float64.(golden["nonlinear_laplacian_u"])
        @test length(nl_lap_u) == length(g_nl)

        @testset "fixture=$(fx["name"])" begin
            for i in 1:N
                scale = max(1.0, abs(nl_lap_u[i]), abs(g_nl[i]))
                @test abs(nl_lap_u[i] - g_nl[i]) <= REL_TOL * scale
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 2D mixed cross-derivative Dxy(u) conformance + 2nd-order convergence
# ---------------------------------------------------------------------------

@testitem "Discretization conformance: rect_2d_mixed_deriv_periodic" setup=[Diffusion2DConformanceHelpers] tags=[:conformance, :discretization] begin
    HARNESS   = joinpath(@__DIR__, "..", "tests", "conformance", "discretization",
                         "rect_2d_mixed_deriv_periodic")
    REPO_ROOT = abspath(joinpath(HARNESS, "..", "..", "..", ".."))
    FIXTURES  = JSON.parsefile(joinpath(HARNESS, "fixtures.json"))
    REL_TOL   = Float64(FIXTURES["tolerance"]["relative"])

    @test get(FIXTURES, "_mol531_sha", nothing) == "35cc9143dc953ac7d3619738bd77b250c1ed162f"

    raw_rule  = JSON.parsefile(joinpath(REPO_ROOT, FIXTURES["rule_path"]))
    rule_name = FIXTURES["rule"]
    rule_body = raw_rule["discretizations"][rule_name]
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw  = rule_lowered["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    for fx in FIXTURES["fixtures"]
        g   = fx["grid"]
        Nx  = Int(g["n_cells_x"])
        Ny  = Int(g["n_cells_y"])
        dx  = (Float64(g["x_end"]) - Float64(g["x_start"])) / Nx
        dy  = (Float64(g["y_end"]) - Float64(g["y_start"])) / Ny
        x0  = Float64(g["x_start"])
        y0  = Float64(g["y_start"])
        xs  = [x0 + (i - 0.5) * dx for i in 1:Nx]
        ys  = [y0 + (j - 0.5) * dy for j in 1:Ny]

        u = [sin(xs[i]) * sin(ys[j]) for i in 1:Nx, j in 1:Ny]

        bindings = Dict{String, Float64}("dx" => dx, "dy" => dy)

        dxy_u = Matrix{Float64}(undef, Nx, Ny)
        for iy in 1:Ny, ix in 1:Nx
            dxy_u[ix, iy] = _eval_replacement_2d(replacement, u, ix, iy, bindings, Nx, Ny)
        end

        golden = JSON.parsefile(joinpath(HARNESS, "golden", "$(fx["name"]).json"))

        @test golden["_mol531_sha"] == "35cc9143dc953ac7d3619738bd77b250c1ed162f"

        g_dxy = golden["dxy_u"]
        @test length(g_dxy) == Nx
        @test length(g_dxy[1]) == Ny

        @testset "fixture=$(fx["name"])" begin
            for ix in 1:Nx, iy in 1:Ny
                v     = dxy_u[ix, iy]
                g_v   = Float64(g_dxy[ix][iy])
                scale = max(1.0, abs(v), abs(g_v))
                @test abs(v - g_v) <= REL_TOL * scale
            end
        end
    end
end

@testitem "Mixed cross-derivative Dxy(u) converges 2nd order on sin(x)*sin(y)" setup=[Diffusion2DConformanceHelpers] tags=[:conformance, :discretization] begin
    # Verify O(dx^2) convergence of the 4-point stencil on u=sin(x)*sin(y) on
    # [0,2π]×[0,2π]. Exact Dxy(u) = cos(x)*cos(y). Grid doubles from N=16→32→64;
    # L∞-error ratio must be ≥ 3.5 (theoretical 4 for exact 2nd order).

    REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
    rule_body = JSON.parsefile(joinpath(REPO_ROOT, "discretizations", "finite_difference",
                               "mixed_deriv_2nd_uniform.json"))["discretizations"]["mixed_deriv_2nd_uniform"]
    rule_lowered = lower_stencil_to_replacement(rule_body)
    repl_raw  = rule_lowered["replacement"]
    replacement = (repl_raw isa AbstractDict && get(repl_raw, "op", nothing) == "arrayop") ?
                  repl_raw["expr"] : repl_raw

    function _linf_error(N)
        dx = 2π / N; dy = 2π / N
        xs = [(i - 0.5) * dx for i in 1:N]
        ys = [(j - 0.5) * dy for j in 1:N]
        u  = [sin(xs[i]) * sin(ys[j]) for i in 1:N, j in 1:N]
        b  = Dict{String,Float64}("dx" => dx, "dy" => dy)
        err = 0.0
        for iy in 1:N, ix in 1:N
            fd  = _eval_replacement_2d(replacement, u, ix, iy, b, N, N)
            ex  = cos(xs[ix]) * cos(ys[iy])
            err = max(err, abs(fd - ex))
        end
        return err
    end

    e16 = _linf_error(16)
    e32 = _linf_error(32)
    e64 = _linf_error(64)
    @test e16 / e32 >= 3.5
    @test e32 / e64 >= 3.5
end
