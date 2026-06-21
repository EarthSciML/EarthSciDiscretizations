# Empirical feasibility probe for declarative spectral operators (esd-6g4.4).
# Tests whether a DENSE reduce arrayop (out[i] = reduce_j body(i,j)) executes
# in the ESS tree-walk backend (the one the ESD Layer-B runners use), and
# whether a Fourier differentiation matrix can be expressed with INLINE AST
# entries (no precomputed matrix => fully declarative).

using EarthSciSerialization
const ESM = EarthSciSerialization

_n(x)   = NumExpr(Float64(x))
_i(x)   = IntExpr(Int64(x))
_v(n)   = VarExpr(String(n))
_op(op, args...; kw...) = OpExpr(String(op), ESM.Expr[args...]; kw...)
_idx(var, idx_exprs...)  = _op("index", _v(var), idx_exprs...)
_D_idx(var, idx_exprs...) = _op("D", _idx(var, idx_exprs...); wrt="t")
_arrayop1d(body, idx, lo, hi) = OpExpr("arrayop", ESM.Expr[];
    output_idx=Any[idx], expr_body=body, ranges=Dict(idx => [lo, hi]))

# arrayop with an output index i AND a reduction index j (j not in output_idx)
_reduce_arrayop(body, i, ilo, ihi, j, jlo, jhi; redop="+") =
    OpExpr("arrayop", ESM.Expr[];
        output_idx=Any[i], expr_body=body,
        ranges=Dict(i => [ilo, ihi], j => [jlo, jhi]), reduce=redop)

function probe(f, name)
    print("[$name] ")
    try
        f()
    catch e
        println("FAIL/THREW: ", sprint(showerror, e))
        return
    end
end

# ----------------------------------------------------------------------
# PROBE 1: dense reduce over a STATE vector — du[i] = Σ_j u[j]
# Tests the contraction mechanism itself in the tree-walk path.
# ----------------------------------------------------------------------
probe("1 dense-sum-reduce") do
    N = 5
    vars = Dict("u" => ModelVariable(StateVariable))
    lhs = _arrayop1d(_D_idx("u", _v("i")), "i", 1, N)
    rhs = _reduce_arrayop(_idx("u", _v("j")), "i", 1, N, "j", 1, N)
    model = ESM.Model(vars, [ESM.Equation(lhs, rhs)])
    ics = Dict("u[$k]" => Float64(k) for k in 1:N)   # u = [1,2,3,4,5], Σ=15
    f!, u0, p, _, vmap = build_evaluator(model; initial_conditions=ics)
    du = similar(u0); f!(du, u0, p, 0.0)
    vals = [du[vmap["u[$k]"]] for k in 1:N]
    ok = all(isapprox.(vals, 15.0; rtol=1e-12))
    println(ok ? "PASS (du=$vals, expect all 15)" : "WRONG (du=$vals, expect all 15)")
end

# ----------------------------------------------------------------------
# PROBE 2: reduction index used as a SCALAR value — du[i] = Σ_j j*u[j]
# Tests whether the reduction index j is bound as a numeric value in
# arithmetic (needed to inline matrix-entry formulas).
# ----------------------------------------------------------------------
probe("2 reduce-index-as-scalar") do
    N = 4
    vars = Dict("u" => ModelVariable(StateVariable))
    lhs = _arrayop1d(_D_idx("u", _v("i")), "i", 1, N)
    body = _op("*", _v("j"), _idx("u", _v("j")))   # j * u[j]
    rhs = _reduce_arrayop(body, "i", 1, N, "j", 1, N)
    model = ESM.Model(vars, [ESM.Equation(lhs, rhs)])
    ics = Dict("u[$k]" => 1.0 for k in 1:N)         # u=1 => Σ j = 1+2+3+4 = 10
    f!, u0, p, _, vmap = build_evaluator(model; initial_conditions=ics)
    du = similar(u0); f!(du, u0, p, 0.0)
    vals = [du[vmap["u[$k]"]] for k in 1:N]
    ok = all(isapprox.(vals, 10.0; rtol=1e-12))
    println(ok ? "PASS (du=$vals, expect all 10)" : "WRONG (du=$vals, expect all 10)")
end

# ----------------------------------------------------------------------
# PROBE 3: FULL Fourier 1st-derivative matrix, inline AST entries.
# Periodic grid x_j = 2π(j-1)/N, N even. u(x)=sin(x) => u'(x)=cos(x).
# D[i,j] = ifelse(i==j, 0, 0.5*(-1)^(i-j)*cot((i-j)π/N)), (-1)^m via cos(πm).
# ----------------------------------------------------------------------
probe("3 fourier-1st-deriv") do
    N = 8
    h = 2π / N
    xs = [(k-1)*h for k in 1:N]
    vars = Dict("u" => ModelVariable(StateVariable))
    lhs = _arrayop1d(_D_idx("u", _v("i")), "i", 1, N)
    m     = _op("-", _v("i"), _v("j"))                       # i - j
    sgn   = _op("cos", _op("*", _n(π), m))                   # (-1)^(i-j)
    theta = _op("*", m, _n(π / N))                           # (i-j)π/N
    cot   = _op("/", _op("cos", theta), _op("sin", theta))   # cot(theta)
    Dij   = _op("*", _n(0.5), sgn, cot)                      # 0.5*(-1)^m*cot
    cond  = _op("==", _v("i"), _v("j"))
    Dsafe = _op("ifelse", cond, _n(0.0), Dij)               # 0 on diagonal
    body  = _op("*", Dsafe, _idx("u", _v("j")))
    rhs   = _reduce_arrayop(body, "i", 1, N, "j", 1, N)
    model = ESM.Model(vars, [ESM.Equation(lhs, rhs)])
    ics = Dict("u[$k]" => sin(xs[k]) for k in 1:N)
    f!, u0, p, _, vmap = build_evaluator(model; initial_conditions=ics)
    du = similar(u0); f!(du, u0, p, 0.0)
    vals = [du[vmap["u[$k]"]] for k in 1:N]
    expect = [cos(xs[k]) for k in 1:N]
    err = maximum(abs.(vals .- expect))
    ok = err < 1e-6
    println(ok ? "PASS (Linf err=$err vs analytic cos)" :
                 "WRONG (Linf err=$err; du=$(round.(vals,digits=4)) expect=$(round.(expect,digits=4)))")
end

println("=== probe complete ===")
