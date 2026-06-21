# Probe 3: (a) Chebyshev 1st-derivative feasibility; (b) Fourier generality at N=16.
using EarthSciSerialization
const ESM = EarthSciSerialization
_n(x) = NumExpr(Float64(x)); _v(n) = VarExpr(String(n))
_op(op, args...; kw...) = OpExpr(String(op), ESM.Expr[args...]; kw...)
_idx(v, is...) = _op("index", _v(v), is...)
_Di(v, is...) = _op("D", _idx(v, is...); wrt = "t")
_ao1(body, idx, lo, hi) = OpExpr("arrayop", ESM.Expr[]; output_idx = Any[idx], expr_body = body, ranges = Dict(idx => [lo, hi]))
_red(body, i, ilo, ihi, j, jlo, jhi) = OpExpr(
    "arrayop", ESM.Expr[]; output_idx = Any[i], expr_body = body,
    ranges = Dict(i => [ilo, ihi], j => [jlo, jhi]), reduce = "+"
)
probe(f, name) = (
    print("[$name] "); try
        f()
    catch e
        println("THREW: ", sprint(showerror, e))
    end
)

# ---- (a) Chebyshev–Gauss–Lobatto 1st derivative, inline AST entries ----
# Points x_k = cos((k-1)π/(N-1)); degree N-1. u=x^3 differentiated exactly => 3x^2.
probe("a chebyshev-1st-deriv") do
    N = 6; Nm1 = N - 1
    xk(k) = cos((k - 1) * π / Nm1)
    xs = [xk(k) for k in 1:N]
    vars = Dict("u" => ModelVariable(StateVariable))
    lhs = _ao1(_Di("u", _v("i")), "i", 1, N)
    # x_i, x_j inline
    xi = _op("cos", _op("*", _op("-", _v("i"), _n(1)), _n(π / Nm1)))
    xj = _op("cos", _op("*", _op("-", _v("j"), _n(1)), _n(π / Nm1)))
    # c_i, c_j (2 at endpoints else 1)
    ci = _op("ifelse", _op("or", _op("==", _v("i"), _n(1)), _op("==", _v("i"), _n(N))), _n(2.0), _n(1.0))
    cj = _op("ifelse", _op("or", _op("==", _v("j"), _n(1)), _op("==", _v("j"), _n(N))), _n(2.0), _n(1.0))
    sgn = _op("cos", _op("*", _n(π), _op("+", _v("i"), _v("j"))))        # (-1)^(i+j)
    offdiag = _op("/", _op("*", _op("/", ci, cj), sgn), _op("-", xi, xj))   # (ci/cj)(-1)^(i+j)/(xi-xj)
    interior = _op("/", _op("-", _n(0.0), xi), _op("*", _n(2.0), _op("-", _n(1.0), _op("*", xi, xi))))  # -xi/(2(1-xi^2))
    corner = _n((2 * Nm1 * Nm1 + 1) / 6)
    diag = _op(
        "ifelse", _op("==", _v("i"), _n(1)), corner,
        _op("ifelse", _op("==", _v("i"), _n(N)), _op("-", _n(0.0), corner), interior)
    )
    Dij = _op("ifelse", _op("==", _v("i"), _v("j")), diag, offdiag)
    body = _op("*", Dij, _idx("u", _v("j")))
    rhs = _red(body, "i", 1, N, "j", 1, N)
    model = ESM.Model(vars, [ESM.Equation(lhs, rhs)])
    ics = Dict("u[$k]" => xs[k]^3 for k in 1:N)
    f!, u0, p, _, vmap = build_evaluator(model; initial_conditions = ics)
    du = similar(u0); f!(du, u0, p, 0.0)
    vals = [du[vmap["u[$k]"]] for k in 1:N]; expect = [3 * xs[k]^2 for k in 1:N]
    err = maximum(abs.(vals .- expect))
    println(
        err < 1.0e-8 ? "PASS (Linf err=$err vs 3x^2, exact poly)" :
            "WRONG (err=$err; du=$(round.(vals, digits = 4)) exp=$(round.(expect, digits = 4)))"
    )
end

# ---- (b) Fourier generality: same matrix at N=16, u=sin(2x) => 2cos(2x) ----
probe("b fourier-N16-wavenumber2") do
    N = 16; h = 2π / N; xs = [(k - 1) * h for k in 1:N]
    vars = Dict("u" => ModelVariable(StateVariable))
    lhs = _ao1(_Di("u", _v("i")), "i", 1, N)
    m = _op("-", _v("i"), _v("j")); sgn = _op("cos", _op("*", _n(π), m))
    th = _op("*", m, _n(π / N)); cot = _op("/", _op("cos", th), _op("sin", th))
    Dij = _op("*", _n(0.5), sgn, cot)
    Ds = _op("ifelse", _op("==", _v("i"), _v("j")), _n(0.0), Dij)
    body = _op("*", Ds, _idx("u", _v("j"))); rhs = _red(body, "i", 1, N, "j", 1, N)
    model = ESM.Model(vars, [ESM.Equation(lhs, rhs)])
    ics = Dict("u[$k]" => sin(2xs[k]) for k in 1:N)
    f!, u0, p, _, vmap = build_evaluator(model; initial_conditions = ics)
    du = similar(u0);f!(du, u0, p, 0.0)
    vals = [du[vmap["u[$k]"]] for k in 1:N];expect = [2cos(2xs[k]) for k in 1:N]
    err = maximum(abs.(vals .- expect))
    println(err < 1.0e-6 ? "PASS (Linf err=$err vs 2cos(2x))" : "WRONG (err=$err)")
end
println("=== probe3 complete ===")
