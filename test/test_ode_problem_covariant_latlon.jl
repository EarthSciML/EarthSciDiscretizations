using Test
using TestItems
import SciMLBase

# ---------------------------------------------------------------------------
# Path A: build_ode_problem routes a covariant-FV (latlon) GDD — one that
# declares the authored declarative discretization rules — through the ESS
# rule-engine pipeline (discretize → build_evaluator), NOT the Path-B PDESystem
# FD chain-rule. This is the G11 production wiring (esd-6g4.10):
#
#   * a latlon GDD WITH a `discretizations` block ⇒ Path A (covariant rules run);
#     a bare-grid latlon GDD ⇒ Path B (unchanged, see test_ode_problem_curvilinear);
#   * the LatLonGrid curvilinear metric is bound as 2-D (nlat, nlon) const_arrays
#     (`_grid_primitive_arrays(::LatLonGrid)`), gathered `index(name, lat, lon)`;
#   * the covariant Laplacian's connection terms gather Jg_xx/Jg_yy/Jg_xe at
#     lat±1 / lon±1 OFFSETS, which resolve through the per-dimension const_array
#     boundary policy (lon :periodic wrap, lat :clamp edge-extend; ess-gj4) —
#     without it build_evaluator would throw E_TREEWALK_CONSTARRAY_OOB at every
#     pole / seam cell.
#
# No imperative covariant operator exists: the discrete Laplace-Beltrami is the
# authored declarative rule, run by the real engine.
#
# ESS stores a [lat,lon] state LAT-FASTEST (Julia column-major of the (nlat,nlon)
# shape): storage index t = lat + (lon-1)*nlat, i.e. lat-row = ((t-1) % nlat)+1.
# ---------------------------------------------------------------------------

@testsnippet CovariantLatLonPathA begin
    import EarthSciDiscretizations as ESD
    import SciMLBase

    _cov_repo() = abspath(joinpath(@__DIR__, ".."))
    _cov_fix(name) = joinpath(_cov_repo(), "test", "fixtures", "curvilinear", name)

    # Build the Path-A ODEProblem for a committed N×N covariant GDD.
    function _cov_build(gdd_name::AbstractString, N::Int)
        esm = _cov_fix("covariant_laplacian_latlon.esm")
        prob, var_map = ESD.build_ode_problem(esm; grid_ref = _cov_fix(gdd_name))
        grid = ESD._latlon(; nlon = N, nlat = N, R = 1.0)
        return prob, var_map, grid
    end

    # State vector holding the degree-1 eigenfunction sin(lat), placed via the
    # var_map's 2-D "u[i,j]" keys (i = lat row, j = lon col) so the test is robust
    # to ESS's internal flat storage order. sin(lat) is lon-independent.
    function _cov_sinlat(var_map, grid, N)
        u = zeros(Float64, N * N)
        for i in 1:N, j in 1:N
            u[var_map["u[$i,$j]"]] = sin(grid.lat_centers[i])
        end
        return u
    end

    # Interior (non-pole lat rows) L∞ error of the discrete operator applied to
    # sin(lat) vs the analytic Laplace-Beltrami −2 sin(lat) (R = 1), evaluated
    # end-to-end through the real ESS engine.
    function _cov_interior_linf(gdd_name::AbstractString, N::Int)
        prob, var_map, grid = _cov_build(gdd_name, N)
        @assert length(prob.u0) == N * N
        u0 = _cov_sinlat(var_map, grid, N)
        du = similar(u0)
        prob.f(du, u0, prob.p, 0.0)
        @assert all(isfinite, du)
        err = 0.0
        for i in 2:(N - 1), j in 1:N      # interior lat rows (exclude poles)
            err = max(err, abs(du[var_map["u[$i,$j]"]] - (-2.0 * sin(grid.lat_centers[i]))))
        end
        return err
    end
end

@testitem "build_ode_problem Path-A: covariant latlon Laplacian builds + runs (ess-gj4 boundary policy resolves pole/seam metric gathers)" setup = [CovariantLatLonPathA] tags = [:covariant_fv] begin
    prob, var_map, _ = _cov_build("covariant_laplacian_latlon.gdd.json", 32)
    @test prob isa SciMLBase.ODEProblem
    @test var_map isa Dict{String, Int}
    @test length(prob.u0) == 32 * 32           # nlat * nlon

    # f! must evaluate over ALL cells — the lat poles (lat±1 metric gathers) and
    # the lon seam (lon±1) — without throwing E_TREEWALK_CONSTARRAY_OOB. A finite
    # RHS here is the end-to-end proof that the LatLonGrid metric is bound as 2-D
    # const_arrays AND the per-dim boundary policy (ess-gj4) is threaded through
    # build_evaluator. (On ESS pre-gj4 this build aborted — esd-6g4.10 verdict.)
    du = similar(prob.u0)
    @test_nowarn prob.f(du, prob.u0, prob.p, 0.0)
    @test all(isfinite, du)
end

@testitem "build_ode_problem Path-A: covariant latlon Laplacian is O(h²) to Laplace-Beltrami on the interior" setup = [CovariantLatLonPathA] tags = [:covariant_fv] begin
    # The interior convergence the rule-conformance test checks via its test-local
    # interpreter, here driven end-to-end through the REAL ESS engine (discretize
    # → build_evaluator + 2-D metric const_arrays + boundary policy). Pole rows
    # carry the documented zero-ghost (state) vs sentinel→self (oracle) divergence
    # (ORACLE_CHARACTERIZATION §6.2) and are excluded from the interior norm.
    e16 = _cov_interior_linf("covariant_laplacian_latlon_n16.gdd.json", 16)
    e32 = _cov_interior_linf("covariant_laplacian_latlon.gdd.json", 32)
    e64 = _cov_interior_linf("covariant_laplacian_latlon_n64.gdd.json", 64)
    @info "covariant latlon Path-A interior L∞" e16 e32 e64 r1 = e16 / e32 r2 = e32 / e64
    @test e16 / e32 >= 3.5          # error quarters per halving ⇒ O(h²)
    @test e32 / e64 >= 3.5
    @test e64 < 2.0e-3
end

# The end-to-end ODE *solve* (a degree-1 eigenfunction decaying as e^{-2t/R²}) is
# exercised as a Layer-C integration fixture, run under the integration test
# target where the test-only ODE solver is on the load path:
#   discretizations/finite_volume/covariant_fv_laplacian_latlon/fixtures/integration/
#   (mms_kind "sin_lat_latlon_laplacian_interior", run via
#   `Pkg.test(; test_args=["integration"])`). Keeping the solver out of this unit
#   file is deliberate — see test/runtests.jl.
