"""
FVCubedSphere discretization: converts a ModelingToolkit PDESystem into an ODEProblem
on a cubed-sphere grid.

Usage:
    sys = PDESystem(eqs, bcs, domains, [t, lon, lat], [u(t, lon, lat)])
    disc = FVCubedSphere(48)
    prob = SciMLBase.discretize(sys, disc)
    sol = solve(prob, Tsit5())
"""

struct FVCubedSphere
    Nc::Int
    Nk::Int
    R::Float64
    Ng::Int
    transport_scheme::Symbol
end

function FVCubedSphere(Nc::Int; Nk::Int=0, R::Float64=6.371e6, Ng::Int=3, transport::Symbol=:upwind)
    FVCubedSphere(Nc, Nk, R, Ng, transport)
end

"""
    SciMLBase.discretize(sys::PDESystem, disc::FVCubedSphere)

Discretize a PDESystem onto a cubed-sphere grid, returning an ODEProblem.

The PDE must have independent variables `(t, x, y)` where the first is time
and the remaining are spatial coordinates on the sphere. Spatial `Differential`
operators are replaced with metric-corrected finite-volume stencil operators.
"""
function SciMLBase.discretize(sys::PDESystem, disc::FVCubedSphere; kwargs...)
    grid = CubedSphereGrid(disc.Nc; R=disc.R, Ng=disc.Ng)
    Nc = grid.Nc

    # Identify time vs spatial IVs
    t_iv = nothing
    spatial_ivs = Num[]
    for iv in sys.ivs
        if identify_dimension(Symbol(iv)) == :t
            t_iv = iv
        else
            push!(spatial_ivs, iv)
        end
    end
    t_iv === nothing && error("No time variable found")

    # Create discrete state array for each dependent variable
    dvs = sys.dvs
    disc_vars = Dict{Any,Any}()
    for dv in dvs
        name = Symbol(Symbolics.tosymbol(dv, escape=false))
        arr = first(Symbolics.@variables $name(mtk_t)[1:6, 1:Nc, 1:Nc])
        disc_vars[dv] = arr
    end

    # Discretize equations
    all_eqs = Symbolics.Equation[]
    for eq in sys.eqs
        new_eqs = discretize_equation(eq, disc_vars, dvs, spatial_ivs, grid, disc)
        append!(all_eqs, new_eqs)
    end

    # Build MTK System
    @named disc_system = System(all_eqs, mtk_t)
    compiled = mtkcompile(disc_system)

    # Project initial conditions
    u0 = build_u0(sys, disc_vars, dvs, grid, compiled)

    # Time span
    tspan = extract_tspan(sys)

    return ODEProblem(compiled, u0, tspan; kwargs...)
end

function extract_tspan(sys::PDESystem)
    for d in sys.domain
        if identify_dimension(Symbol(d.variables)) == :t
            return (Float64(d.domain.left), Float64(d.domain.right))
        end
    end
    error("No time domain found")
end

function build_u0(sys::PDESystem, disc_vars, dvs, grid, compiled)
    Nc = grid.Nc
    u0 = Pair[]

    spatial_ivs = [iv for iv in sys.ivs if identify_dimension(Symbol(iv)) != :t]
    tspan = extract_tspan(sys)
    t0 = tspan[1]

    for dv in dvs
        arr = disc_vars[dv]
        ic_found = false

        for bc in sys.bcs
            if _is_initial_condition(bc, dv, t0)
                rhs = bc.rhs
                for p in 1:6, i in 1:Nc, j in 1:Nc
                    val = _eval_ic(rhs, spatial_ivs, grid, p, i, j)
                    push!(u0, arr[p, i, j] => val)
                end
                ic_found = true
                break
            end
        end

        if !ic_found
            # Default to zero if no IC found for this DV
            for p in 1:6, i in 1:Nc, j in 1:Nc
                push!(u0, arr[p, i, j] => 0.0)
            end
        end
    end

    return u0
end

"""
Check if a boundary condition is an initial condition for a given DV.
An IC has the form `dv(t0, x, y) ~ expr(x, y)` where the first argument
of the DV call equals the start of the time domain.
"""
function _is_initial_condition(bc, dv, t0)
    lhs = unwrap(bc.lhs)
    if !iscall(lhs)
        return false
    end
    args = arguments(lhs)
    # Check name match
    lhs_name = Symbol(Symbolics.tosymbol(wrap(lhs), escape=false))
    dv_name = Symbol(Symbolics.tosymbol(dv, escape=false))
    if lhs_name != dv_name
        return false
    end
    # Check first argument is t0
    if length(args) >= 1
        t_arg = args[1]
        t_val = Symbolics.value(t_arg)
        if t_val isa Number && isapprox(Float64(t_val), t0)
            return true
        end
    end
    return false
end

function _eval_ic(rhs, spatial_ivs, grid, p, i, j)
    try
        subs = Dict{Any,Any}()
        if length(spatial_ivs) >= 1
            subs[spatial_ivs[1]] = grid.lon[p, i, j]
        end
        if length(spatial_ivs) >= 2
            subs[spatial_ivs[2]] = grid.lat[p, i, j]
        end
        v = Symbolics.value(Symbolics.substitute(rhs, subs))
        return v isa Number ? Float64(v) : Float64(eval(Symbolics.toexpr(v)))
    catch
        return 0.0
    end
end
