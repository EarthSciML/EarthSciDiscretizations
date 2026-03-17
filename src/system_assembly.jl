"""
System assembly helpers.
"""

function extract_tspan(domains)
    for d in domains
        hasproperty(d, :variables) && Symbol(d.variables) == :t &&
            return (Float64(d.domain.left), Float64(d.domain.right))
    end
    error("No time domain found")
end
