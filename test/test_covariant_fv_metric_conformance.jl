# ---------------------------------------------------------------------------
# Covariant FV metric operator — oracle conformance (esd-zk9.2 / R2)
#
# The declarative arrayop-einsum rules
#   discretizations/finite_volume/covariant_fv_laplacian_latlon.json
#   discretizations/finite_volume/covariant_fv_gradient_latlon.json
# must be tolerance-identical to the R1 oracle goldens
#   discretizations/finite_volume/covariant_fv_metric/fixtures/oracle/*.json
# (frozen Float64 outputs of the ESS grid_assembly.jl oracle, esd-zk9.1).
#
# This test parses the AUTHORED rule JSON and evaluates its `replacement`
# expression per cell through a small test-local interpreter that resolves the
# `index(...)` gathers under the oracle's exact boundary policy (lon periodic;
# lat poles use the oracle's sentinel→self / reuse-center convention,
# ORACLE_CHARACTERIZATION.md §1.3/§6.2) and delegates every index-free
# arithmetic subtree to the REAL ESS evaluator via `eval_coeff`. So the metric
# algebra (g^{ij} contraction, the J·g^{ij} connection-term centered
# differences, the corner cross-derivative) is exercised through the engine's
# arithmetic, not a shadow reimplementation — only the gather indexing and the
# linear-combination structure are test-local, mirroring the established
# test_diffusion_discretization_conformance.jl pattern.
# ---------------------------------------------------------------------------

@testsnippet CovariantFVConformanceHelpers begin
    using EarthSciDiscretizations: eval_coeff
    using JSON

    _has_index_cov(node)::Bool = begin
        node isa AbstractDict || return false
        get(node, "op", nothing) == "index" && return true
        any(_has_index_cov(a) for a in get(node, "args", Any[]))
    end

    # An index-axis subscript is either the bare axis name (offset 0) or
    # {"op":"+"/"-","args":[axis, k]} with integer k.
    function _axis_offset_cov(expr, axis::AbstractString)::Int
        if expr isa AbstractString
            String(expr) == axis || error("axis name mismatch: '$expr' vs '$axis'")
            return 0
        end
        expr isa AbstractDict || error("unexpected axis subscript $(typeof(expr))")
        a = expr["args"]
        String(a[1]) == axis || error("axis lhs '$(a[1])' != '$axis'")
        k = Int(a[2])
        op = String(expr["op"])
        op == "+" && return k
        op == "-" && return -k
        error("unsupported axis op '$op'")
    end

    # Resolve a stencil neighbour of cell (i,j) [1-based: i along ξ=lon ∈ 1:Nx,
    # j along η=lat ∈ 1:Ny] offset by (di,dj) to a 1-based flat index
    # c = i + (j-1)*Nx. Periodic axes wrap; a non-periodic axis stepped
    # off-domain maps the WHOLE neighbour to self (sentinel→self), exactly as
    # the grid_assembly oracle does at the lat poles.
    function _nb_flat_cov(
            i::Int, j::Int, di::Int, dj::Int,
            Nx::Int, Ny::Int, pxi::Bool, peta::Bool
        )::Int
        self = i + (j - 1) * Nx
        jp = j + dj
        if jp < 1 || jp > Ny
            peta ? (jp = mod(jp - 1, Ny) + 1) : return self
        end
        ip = i + di
        if ip < 1 || ip > Nx
            pxi ? (ip = mod(ip - 1, Nx) + 1) : return self
        end
        return ip + (jp - 1) * Nx
    end

    # Evaluate a rule `replacement` node at cell (i,j). `arrays` maps each field
    # / metric symbol ("$u", "g_xx", "Jg_xe", "dxi_dt1", …) to its flat per-cell
    # Vector; `bindings` carries the scalar params (dlon, dlat).
    function _eval_cov(
            node, arrays::Dict{String, Vector{Float64}},
            i::Int, j::Int, Nx::Int, Ny::Int, pxi::Bool, peta::Bool,
            bindings::Dict{String, Float64}
        )::Float64
        node isa Number && return Float64(node)
        if node isa AbstractString
            s = String(node)
            haskey(bindings, s) && return bindings[s]
            error("unresolved variable '$s'")
        end
        node isa AbstractDict || error("unexpected node type $(typeof(node))")
        op = String(node["op"])
        args = node["args"]
        if op == "index"
            name = String(args[1])
            dj = _axis_offset_cov(args[2], "lat")   # η / lat subscript
            di = _axis_offset_cov(args[3], "lon")   # ξ / lon subscript
            flat = _nb_flat_cov(i, j, di, dj, Nx, Ny, pxi, peta)
            return arrays[name][flat]
        end
        # Index-free subtree → real ESS arithmetic (DomainError/Unbound etc. all
        # propagate from the engine, not a local reimplementation).
        _has_index_cov(node) || return eval_coeff(node, bindings)
        ev(a) = _eval_cov(a, arrays, i, j, Nx, Ny, pxi, peta, bindings)
        op == "+" && return length(args) == 1 ? ev(args[1]) : sum(ev(a) for a in args)
        op == "-" && return length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2])
        op == "*" && return length(args) == 1 ? ev(args[1]) : prod(ev(a) for a in args)
        op == "/" && return ev(args[1]) / ev(args[2])
        error("unsupported op '$op'")
    end

    function _load_repl_cov(repo_root, file, name)
        body = JSON.parsefile(
            joinpath(
                repo_root, "discretizations", "finite_volume", file
            )
        )["discretizations"][name]
        r = body["replacement"]
        return (r isa AbstractDict && get(r, "op", nothing) == "arrayop") ? r["expr"] : r
    end
end

@testitem "Covariant FV metric: declarative rule is tolerance-identical to grid_assembly oracle goldens (esd-zk9.2)" setup = [CovariantFVConformanceHelpers] tags = [:conformance, :covariant_fv] begin
    using JSON

    REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
    ORACLE_DIR = joinpath(
        REPO_ROOT, "discretizations", "finite_volume",
        "covariant_fv_metric", "fixtures", "oracle"
    )

    lap_repl = _load_repl_cov(REPO_ROOT, "covariant_fv_laplacian_latlon.json", "covariant_fv_laplacian_latlon")
    t1_repl = _load_repl_cov(REPO_ROOT, "covariant_fv_gradient_latlon.json", "covariant_fv_gradient_latlon_t1")
    t2_repl = _load_repl_cov(REPO_ROOT, "covariant_fv_gradient_latlon.json", "covariant_fv_gradient_latlon_t2")

    # Oracle and rule share Float64 arithmetic; the only residual is operation-
    # ordering roundoff (≈1e-15 absolute on |result| up to ~1.7). A relative-or-
    # absolute window 10 orders below the signal is "tolerance-identical".
    TOL = 1.0e-10

    for case in ("latlon_small", "nonorthogonal_small")
        d = JSON.parsefile(joinpath(ORACLE_DIR, "$case.json"))
        g = d["grid"]
        Nx = Int(g["Nx"]); Ny = Int(g["Ny"]); nc = Nx * Ny
        dlon = Float64(g["dxi"]); dlat = Float64(g["deta"])
        pxi = Bool(g["periodic_xi"]); peta = Bool(g["periodic_eta"])

        m = d["metric"]
        gxx = Float64.(m["ginv_xixi"]); gyy = Float64.(m["ginv_etaeta"])
        gxe = Float64.(m["ginv_xieta"]); J = Float64.(m["jacobian_J"])
        phi = Float64.(d["field_phi"])
        cj = d["coord_jacobian"]

        arrays = Dict{String, Vector{Float64}}(
            "\$u" => phi,
            "g_xx" => gxx, "g_yy" => gyy, "g_xe" => gxe,
            "invJ" => Float64[1.0 / J[k] for k in 1:nc],
            "Jg_xx" => Float64[J[k] * gxx[k] for k in 1:nc],
            "Jg_yy" => Float64[J[k] * gyy[k] for k in 1:nc],
            "Jg_xe" => Float64[J[k] * gxe[k] for k in 1:nc],
            "dxi_dt1" => Float64.(cj["dxi_dt1"]), "deta_dt1" => Float64.(cj["deta_dt1"]),
            "dxi_dt2" => Float64.(cj["dxi_dt2"]), "deta_dt2" => Float64.(cj["deta_dt2"]),
        )
        bindings = Dict{String, Float64}("dlon" => dlon, "dlat" => dlat)

        lap_o = Float64.(d["laplacian"]["applied_result"])
        t1_o = Float64.(d["gradient"]["applied_result_t1"])
        t2_o = Float64.(d["gradient"]["applied_result_t2"])

        lap = Vector{Float64}(undef, nc)
        t1 = Vector{Float64}(undef, nc)
        t2 = Vector{Float64}(undef, nc)

        @testset "case=$case (Nx=$Nx, Ny=$Ny, periodic_eta=$peta)" begin
            for j in 1:Ny, i in 1:Nx
                c = i + (j - 1) * Nx
                lap[c] = _eval_cov(lap_repl, arrays, i, j, Nx, Ny, pxi, peta, bindings)
                t1[c] = _eval_cov(t1_repl, arrays, i, j, Nx, Ny, pxi, peta, bindings)
                t2[c] = _eval_cov(t2_repl, arrays, i, j, Nx, Ny, pxi, peta, bindings)
                @test abs(lap[c] - lap_o[c]) <= TOL * max(1.0, abs(lap[c]), abs(lap_o[c]))
                @test abs(t1[c] - t1_o[c]) <= TOL * max(1.0, abs(t1[c]), abs(t1_o[c]))
                @test abs(t2[c] - t2_o[c]) <= TOL * max(1.0, abs(t2[c]), abs(t2_o[c]))
            end

            # Whole-operator checksums (ORACLE_CHARACTERIZATION fixtures schema).
            ck = d["checksums"]
            @test abs(sum(lap) - Float64(ck["sum_du_lap"])) <= 1.0e-9
            @test abs(sum(abs, lap) - Float64(ck["sumabs_du_lap"])) <= TOL * abs(Float64(ck["sumabs_du_lap"]))
            @test abs(maximum(abs, lap) - Float64(ck["maxabs_du_lap"])) <= TOL * abs(Float64(ck["maxabs_du_lap"]))
            @test abs(sum(abs, t1) - Float64(ck["sumabs_grad_t1"])) <= TOL * abs(Float64(ck["sumabs_grad_t1"]))
            @test abs(sum(abs, t2) - Float64(ck["sumabs_grad_t2"])) <= TOL * abs(Float64(ck["sumabs_grad_t2"]))
        end
    end
end

@testitem "Covariant FV metric: PDE solve/verify — rule + LatLonGrid binding hook converges to Laplace-Beltrami O(h^2) (esd-zk9.2)" setup = [CovariantFVConformanceHelpers] tags = [:conformance, :covariant_fv] begin
    using EarthSciDiscretizations: n_cells, cell_centers
    import EarthSciDiscretizations as ESD
    using JSON

    # End-to-end PDE solve/verify driven by the declarative rule: build the
    # metric const_arrays exactly as production would — the host-side binding
    # hook `_grid_primitive_arrays(::LatLonGrid)` (esd-zk9.2 §6.1, MPAS/DUO
    # precedent) — feed them to the authored covariant Laplacian replacement,
    # and verify the discrete operator applied to the degree-1 spherical
    # eigenfunction f = sin(lat) converges at O(h²) to the analytic
    # Laplace-Beltrami result ∇²f = −2 sin(lat)/R². The lat-pole rows carry the
    # documented sentinel→self vs zero-ghost boundary divergence (§6.2) and are
    # excluded from the interior convergence norm.
    REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
    lap_repl = _load_repl_cov(REPO_ROOT, "covariant_fv_laplacian_latlon.json", "covariant_fv_laplacian_latlon")
    R = 1.0

    function interior_linf(nlon::Int, nlat::Int)
        grid = ESD._latlon(; nlon = nlon, nlat = nlat, R = R)
        nc = n_cells(grid)
        # The binding hook now returns (nlat, nlon) metric MATRICES (esd-6g4.10):
        # ESS gathers index(name, lat, lon) as a 2-D const_array. This test-local
        # interpreter resolves gathers against a lon-fastest FLAT layout
        # (c = lon + (lat-1)*nlon), so flatten each matrix back to that order
        # (inverse of the production reshape permutedims(reshape(v,nlon,nlat),(2,1))).
        prim = ESD._grid_primitive_arrays(grid)   # the binding hook under test
        _flat_lonfastest(M) = vec(permutedims(M, (2, 1)))
        cc = cell_centers(grid)
        lat = cc.lat
        lats = sort(unique(lat))
        dlon = 2π / nlon
        dlat = lats[2] - lats[1]
        # regular variant ⇒ uniform lat rows; a non-uniform grid would void the
        # single-scalar dlat and is caught here rather than as a silent miss.
        @test all(abs((lats[k + 1] - lats[k]) - dlat) < 1.0e-10 for k in 1:(length(lats) - 1))

        arrays = Dict{String, Vector{Float64}}("\$u" => Float64[sin(lat[c]) for c in 1:nc])
        for k in ("g_xx", "g_yy", "g_xe", "invJ", "Jg_xx", "Jg_yy", "Jg_xe")
            arrays[k] = _flat_lonfastest(Float64.(prim[k]))
        end
        bindings = Dict{String, Float64}("dlon" => dlon, "dlat" => dlat)

        err = 0.0
        for j in 2:(nlat - 1), i in 1:nlon   # interior lat rows (exclude poles)
            c = i + (j - 1) * nlon
            lap = _eval_cov(lap_repl, arrays, i, j, nlon, nlat, true, false, bindings)
            exact = -2.0 * sin(lat[c]) / (R * R)
            err = max(err, abs(lap - exact))
        end
        return err
    end

    e16 = interior_linf(16, 16)
    e32 = interior_linf(32, 32)
    e64 = interior_linf(64, 64)
    @test e16 / e32 >= 3.5          # O(h²) ⇒ error quarters per halving
    @test e32 / e64 >= 3.5
    @test e64 < 2.0e-3
end
