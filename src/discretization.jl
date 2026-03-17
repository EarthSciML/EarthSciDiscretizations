"""
FVCubedSphere discretization specification.
"""

struct FVCubedSphere
    Nc::Int; Nk::Int; R::Float64; Ng::Int; transport_scheme::Symbol
end

function FVCubedSphere(Nc::Int; Nk::Int = 0, R::Float64 = 6.371e6, Ng::Int = 3, transport::Symbol = :upwind)
    FVCubedSphere(Nc, Nk, R, Ng, transport)
end

function infer_var_locations(dvs)
    [Symbol(dv) => CellCenter for dv in dvs]
end
