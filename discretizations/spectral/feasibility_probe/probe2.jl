# Probe 2: does a Fourier reduce-arrayop RULE apply through the real
# EarthSciSerialization.discretize pipeline, and can the reduction range be
# tied to the (dynamic) grid dimension size so the rule is GENERAL?

using EarthSciSerialization
import JSON
const ESM = EarthSciSerialization

# Build the Fourier 1st-derivative matrix entry D[i,j] as a JSON-AST dict.
# D[i,j] = ifelse(i==j, 0, 0.5*(-1)^(i-j)*cot((i-j)*pi/N)); cot=cos/sin;
# (-1)^m = cos(pi*m). pi_over_N is passed as a binding so the rule stays
# grid-agnostic; we also try a literal version.
function Dij_ast(pi_over_N_expr)
    m     = Dict("op"=>"-", "args"=>["i","j"])
    sgn   = Dict("op"=>"cos", "args"=>[Dict("op"=>"*", "args"=>[Dict("op"=>"pi","args"=>[]), m])])
    theta = Dict("op"=>"*", "args"=>[m, pi_over_N_expr])
    cot   = Dict("op"=>"/", "args"=>[Dict("op"=>"cos","args"=>[theta]), Dict("op"=>"sin","args"=>[theta])])
    Dij   = Dict("op"=>"*", "args"=>[0.5, sgn, cot])
    cond  = Dict("op"=>"==", "args"=>["i","j"])
    Dict("op"=>"ifelse", "args"=>[cond, 0.0, Dij])
end

function build_input(N::Int; reduce_hi)
    body = Dict("op"=>"*", "args"=>[
        Dij_ast(Dict("op"=>"/", "args"=>[Dict("op"=>"pi","args"=>[]), "Nx"])),
        Dict("op"=>"index", "args"=>["\$u", "j"]),
    ])
    replacement = Dict(
        "op"=>"arrayop", "reduce"=>"+",
        "ranges"=>Dict("j"=>[1, reduce_hi]),
        "expr"=>body,
    )
    Dict(
        "esm"=>"0.4.0",
        "grids"=>Dict("g"=>Dict(
            "family"=>"cartesian",
            "dimensions"=>[Dict("name"=>"x","periodic"=>true,"size"=>N,"spacing"=>"uniform")],
        )),
        "metadata"=>Dict("name"=>"fourier_diff_1d_probe"),
        "models"=>Dict("M"=>Dict(
            "grid"=>"g",
            "variables"=>Dict(
                "u"=>Dict("type"=>"state","location"=>"cell_center","shape"=>["x"],"default"=>0.0,"units"=>"1"),
                "Nx"=>Dict("type"=>"parameter","default"=>Float64(N),"units"=>"1"),
            ),
            "equations"=>[Dict(
                "lhs"=>Dict("op"=>"D","wrt"=>"t","args"=>["u"]),
                "rhs"=>Dict("op"=>"grad","args"=>["u"],"dim"=>"x"),
            )],
        )),
        "rules"=>[Dict(
            "name"=>"fourier_diff_1d",
            "pattern"=>Dict("op"=>"grad","args"=>["\$u"],"dim"=>"x"),
            "replacement"=>replacement,
        )],
    )
end

function try_case(label, N, reduce_hi)
    print("[$label] ")
    input = build_input(N; reduce_hi=reduce_hi)
    out = try
        ESM.discretize(input; lift_1d_arrayop=true)
    catch e
        println("discretize THREW: ", sprint(showerror, e)); return nothing
    end
    # Did the grad get rewritten to an arrayop with a reduce?
    eqs = out["models"]["M"]["equations"]
    rhs = eqs[1]["rhs"]
    println("discretize OK; rhs op=", get(rhs,"op","?"),
            " reduce=", get(rhs,"reduce","<none>"),
            " ranges=", get(rhs,"ranges","<none>"))
    return out
end

println("--- Case A: reduction range hard-coded to grid size (N=8) ---")
outA = try_case("A hardcoded [1,8]", 8, 8)

println("--- Case B: reduction range as symbolic grid-size string \"Nx\" ---")
outB = try_case("B symbolic \"Nx\"", 8, "Nx")

println("--- Case C: reduction range as dimension name \"x\" ---")
outC = try_case("C dim-name \"x\"", 8, "x")

# If Case A discretized, evaluate it to confirm correctness end-to-end.
if outA !== nothing
    print("[A eval] ")
    try
        N = 8; h = 2π/N; xs = [(k-1)*h for k in 1:N]
        f!, u0, p, _, vmap = ESM.build_evaluator(outA;
            initial_conditions=Dict("u[$k]"=>sin(xs[k]) for k in 1:N),
            const_arrays=Dict{String,AbstractArray{Float64}}())
        du = similar(u0); f!(du, u0, p, 0.0)
        vals = [du[vmap["u[$k]"]] for k in 1:N]
        expect = [cos(xs[k]) for k in 1:N]
        err = maximum(abs.(vals .- expect))
        println(err < 1e-6 ? "PASS through real pipeline (Linf err=$err)" : "WRONG (err=$err)")
    catch e
        println("eval THREW: ", sprint(showerror, e))
    end
end
println("=== probe2 complete ===")
