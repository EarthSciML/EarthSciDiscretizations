# Verify the GENERAL Fourier rule (reduction range tied to the grid, not
# hardcoded) applies + evaluates to machine precision at multiple grid sizes.
using EarthSciSerialization
import JSON
const ESM = EarthSciSerialization

function Dij_ast()
    m     = Dict("op"=>"-", "args"=>["i","j"])
    sgn   = Dict("op"=>"cos", "args"=>[Dict("op"=>"*", "args"=>[Dict("op"=>"pi","args"=>[]), m])])
    theta = Dict("op"=>"*", "args"=>[m, Dict("op"=>"/", "args"=>[Dict("op"=>"pi","args"=>[]), "Nx"])])
    cot   = Dict("op"=>"/", "args"=>[Dict("op"=>"cos","args"=>[theta]), Dict("op"=>"sin","args"=>[theta])])
    Dij   = Dict("op"=>"*", "args"=>[0.5, sgn, cot])
    cond  = Dict("op"=>"==", "args"=>["i","j"])
    Dict("op"=>"ifelse", "args"=>[cond, 0.0, Dij])
end

function build_input(N::Int, reduce_hi)
    body = Dict("op"=>"*", "args"=>[Dij_ast(), Dict("op"=>"index", "args"=>["\$u","j"])])
    repl = Dict("op"=>"arrayop","reduce"=>"+","ranges"=>Dict("j"=>[1,reduce_hi]),"expr"=>body)
    Dict("esm"=>"0.4.0",
      "grids"=>Dict("g"=>Dict("family"=>"cartesian",
        "dimensions"=>[Dict("name"=>"x","periodic"=>true,"size"=>N,"spacing"=>"uniform")])),
      "metadata"=>Dict("name"=>"fourier_diff_1d"),
      "models"=>Dict("M"=>Dict("grid"=>"g",
        "variables"=>Dict(
          "u"=>Dict("type"=>"state","location"=>"cell_center","shape"=>["x"],"default"=>0.0,"units"=>"1"),
          "Nx"=>Dict("type"=>"parameter","default"=>Float64(N),"units"=>"1")),
        "equations"=>[Dict("lhs"=>Dict("op"=>"D","wrt"=>"t","args"=>["u"]),
                           "rhs"=>Dict("op"=>"grad","args"=>["u"],"dim"=>"x"))])),
      "rules"=>[Dict("name"=>"fourier_diff_1d",
        "pattern"=>Dict("op"=>"grad","args"=>["\$u"],"dim"=>"x"),
        "replacement"=>repl)])
end

function run(N, reduce_hi, label)
    print("[$label N=$N] ")
    input = build_input(N, reduce_hi)
    out = try ESM.discretize(input; lift_1d_arrayop=true)
          catch e; println("discretize THREW: ", sprint(showerror,e)); return end
    h=2π/N; xs=[(k-1)*h for k in 1:N]
    try
        f!,u0,p,_,vmap = ESM.build_evaluator(out;
            initial_conditions=Dict("u[$k]"=>sin(xs[k]) for k in 1:N),
            const_arrays=Dict{String,AbstractArray{Float64}}())
        du=similar(u0); f!(du,u0,p,0.0)
        vals=[du[vmap["u[$k]"]] for k in 1:N]; expect=[cos(xs[k]) for k in 1:N]
        err=maximum(abs.(vals.-expect))
        println(err<1e-6 ? "PASS (Linf err=$err)" : "WRONG (err=$err)")
    catch e; println("eval THREW: ", sprint(showerror,e)) end
end

# General rule: reduction range upper bound = "Nx" (grid-size parameter)
run(8,  "Nx", "param-Nx")
run(12, "Nx", "param-Nx")
run(16, "Nx", "param-Nx")
# Alternative: dimension-name token "x"
run(8,  "x",  "dim-x")
run(12, "x",  "dim-x")
println("=== verify_general complete ===")
