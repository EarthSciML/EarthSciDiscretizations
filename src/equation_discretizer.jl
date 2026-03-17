"""
Dimension identification helpers for the discretization pipeline.
"""

function identify_dimension(iv)
    name = Symbol(iv)
    name == :t && return :t
    name in (:lon, :λ, :x, :ξ, :xi) && return :xi
    name in (:lat, :φ, :y, :η, :eta) && return :eta
    name in (:z, :p, :σ, :k) && return :vertical
    error("Cannot identify grid dimension for variable: $name")
end
