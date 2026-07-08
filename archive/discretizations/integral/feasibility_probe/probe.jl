# Integral / PIDE declarative-feasibility probe (esd-6g4.14, "G15").
#
# Drives the OFFICIAL ESS tree-walk evaluator (`discretize` / `build_evaluator`)
# — the same canonical pipeline the ESD walker uses — to establish that a domain
# integral CAN be expressed as a declarative `reduce` (sum_product) arrayop RULE.
#
# Verdict: FEASIBLE. The shipped rule is `../midpoint_1d.json`; see
# ../INTEGRAL_FEASIBILITY.md. Six checks:
#   1. The ESS `integral` OP is already grid-agnostic (existing host quadrature).
#   2. A scalar `reduce` arrayop reproduces the midpoint integral at a CONCRETE N.
#   3. The NAIVE grid-tied bound (param "Nx" / dim "x") throws
#      E_TREEWALK_UNBOUND_LOOP_VAR — the apparent GAP A blocker.
#   4. *** KEY: a const_array-supplied bound `index(size_x,1)` BYPASSES GAP A ***
#      — the reduction tracks the grid size as DATA, grid-agnostic, err=0.0.
#   5. The full RULE path: the shipped midpoint_1d pattern/replacement rewrites
#      the `integral` op → arrayop via discretize and evaluates correctly.
#   6. The compact-support operator shape (per-cell bound) — the integral
#      OPERATOR generalisation, also feasible in shape.
using EarthSciSerialization
import JSON
const ESM = EarthSciSerialization
idx(a, ks...) = Dict("op" => "index", "args" => Any[a, ks...])

function wrap_esm(N::Int, rhs; extra_vars = Dict{String, Any}(), rules = nothing)
    vars = Dict{String, Any}(
        "u" => Dict("type" => "state", "location" => "cell_center", "shape" => ["x"], "default" => 0.0, "units" => "1")
    )
    merge!(vars, extra_vars)
    eq = Dict("lhs" => Dict("op" => "D", "wrt" => "t", "args" => ["u"]), "rhs" => rhs)
    d = Dict(
        "esm" => "0.4.0", "metadata" => Dict("name" => "integral_probe"),
        "grids" => Dict(
            "g" => Dict(
                "family" => "cartesian",
                "dimensions" => [Dict("name" => "x", "periodic" => true, "size" => N, "spacing" => "uniform")]
            )
        ),
        "models" => Dict("M" => Dict("grid" => "g", "variables" => vars, "equations" => [eq]))
    )
    rules !== nothing && (d["rules"] = rules)
    return d
end

function check(label, input, N, expect; const_arrays = Dict{String, AbstractArray{Float64}}())
    print("[$label] ")
    out = try
        ESM.discretize(input; lift_1d_arrayop = true)
    catch e
        println("discretize THREW: ", first(split(sprint(showerror, e), '\n'))); return
    end
    return try
        f!, u0, p, _, vmap = ESM.build_evaluator(
            out;
            initial_conditions = Dict("u[$k]" => 1.0 for k in 1:N), const_arrays = const_arrays
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        vals = [du[vmap["u[$k]"]] for k in 1:N]; err = maximum(abs.(vals .- expect))
        println(err < 1.0e-10 ? "PASS  du≡$(expect)  (Linf err=$err)" : "WRONG (err=$err, vals[1]=$(vals[1]))")
    catch e
        println("eval THREW: ", first(split(sprint(showerror, e), '\n')))
    end
end

# --- Check 1: ESS `integral` op (existing host quadrature, grid-agnostic) ---
function integral_op(N)
    integral = Dict(
        "op" => "integral", "args" => ["u"], "var" => "x",
        "lower" => Dict("op" => "const", "value" => 0.0), "upper" => Dict("op" => "const", "value" => 1.0)
    )
    return wrap_esm(
        N, Dict("op" => "neg", "args" => [integral]);
        extra_vars = Dict("dx" => Dict("type" => "parameter", "default" => 1.0 / N, "units" => "1"))
    )
end

# --- Checks 2 & 3: scalar reduce-arrayop integral, du/dt = -Σ_j (h·u[j]) ---
function arrayop_integral(N, reduce_hi)
    body = Dict("op" => "*", "args" => [1.0 / N, idx("u", "j")])
    repl = Dict("op" => "arrayop", "reduce" => "+", "ranges" => Dict("j" => [1, reduce_hi]), "expr" => body)
    return wrap_esm(
        N, Dict("op" => "neg", "args" => [repl]);
        extra_vars = Dict("Nx" => Dict("type" => "parameter", "default" => Float64(N), "units" => "1"))
    )
end

# --- Check 4: const_array-supplied bound index(size_x,1) — GAP A BYPASS ---
function const_sized_integral(N)
    body = Dict("op" => "*", "args" => [idx("w", "j"), idx("u", "j")])
    repl = Dict(
        "op" => "arrayop", "reduce" => "+", "ranges" => Dict("j" => [1, idx("size_x", 1)]),
        "expr" => body, "args" => ["w", "size_x"]
    )
    return wrap_esm(N, Dict("op" => "neg", "args" => [repl]))
end

# --- Check 5: full RULE path — the shipped midpoint_1d rule rewrites integral op ---
function rule_path_integral(N)
    integral = Dict(
        "op" => "integral", "args" => ["u"], "var" => "x",
        "lower" => Dict("op" => "const", "value" => 0.0), "upper" => Dict("op" => "const", "value" => 1.0)
    )
    repl = Dict(
        "op" => "*", "args" => [
            "dx",
            Dict(
                "op" => "arrayop", "reduce" => "+", "ranges" => Dict("j" => [1, idx("size_x", 1)]),
                "expr" => idx("u", "j"), "args" => ["size_x"]
            ),
        ]
    )
    rule = Dict(
        "name" => "midpoint_1d",
        "pattern" => Dict(
            "op" => "integral", "args" => ["\$u"], "var" => "x",
            "lower" => Dict("op" => "const", "value" => 0.0), "upper" => Dict("op" => "const", "value" => 1.0)
        ),
        "replacement" => repl
    )
    return wrap_esm(
        N, Dict("op" => "neg", "args" => [integral]);
        extra_vars = Dict("dx" => Dict("type" => "parameter", "default" => 1.0 / N, "units" => "1")), rules = [rule]
    )
end

# --- Check 6: compact-support operator (nn_diffusion shape, per-cell bound) ---
function compact_integral(N)
    body = Dict("op" => "*", "args" => [idx("w", "i", "k"), idx("u", idx("supp", "i", "k"))])
    repl = Dict(
        "op" => "arrayop", "reduce" => "+",
        "ranges" => Dict("k" => [1, idx("n_supp", "i")]), "expr" => body, "args" => ["w", "supp", "n_supp"]
    )
    return wrap_esm(N, Dict("op" => "neg", "args" => [repl]))
end

println("=== Integral/PIDE declarative-feasibility probe (esd-6g4.14) — verdict: FEASIBLE ===")
println("--- Check 1: ESS `integral` op grid-agnostic (existing host quadrature) ---")
for N in (8, 16, 32)
    check("1: integral OP      N=$N", integral_op(N), N, -1.0)
end
println("--- Check 2: scalar reduce-arrayop integral, CONCRETE bound (mechanism) ---")
for N in (8, 16)
    check("2: arrayop concrete  N=$N", arrayop_integral(N, N), N, -1.0)
end
println("--- Check 3: NAIVE grid-tied symbolic bound fails (the apparent GAP A) ---")
check("3: arrayop param-Nx  N=8 ", arrayop_integral(8, "Nx"), 8, -1.0)
check("3: arrayop dim-x     N=16", arrayop_integral(16, "x"), 16, -1.0)
println("--- Check 4: *** const_array bound index(size_x,1) BYPASSES GAP A *** ---")
for N in (8, 16, 32)
    check(
        "4: const-array bound N=$N", const_sized_integral(N), N, -1.0;
        const_arrays = Dict{String, AbstractArray{Float64}}("w" => fill(1.0 / N, N), "size_x" => Float64[N])
    )
end
println("--- Check 5: full RULE path (shipped midpoint_1d rewrites the integral op) ---")
for N in (8, 16, 32)
    check(
        "5: rule path        N=$N", rule_path_integral(N), N, -1.0;
        const_arrays = Dict{String, AbstractArray{Float64}}("size_x" => Float64[N])
    )
end
println("--- Check 6: compact-support operator shape (integral OPERATOR generalisation) ---")
for N in (8, 16)
    h = 1.0 / N; supp = Array{Float64}(undef, N, 3); w = fill(h, N, 3); nsupp = fill(3.0, N)
    for i in 1:N
        supp[i, 1] = mod1(i - 1, N); supp[i, 2] = i; supp[i, 3] = mod1(i + 1, N)
    end
    check(
        "6: compact OP       N=$N", compact_integral(N), N, -3.0 * h;
        const_arrays = Dict{String, AbstractArray{Float64}}("supp" => supp, "w" => w, "n_supp" => nsupp)
    )
end
println("=== probe complete ===")
