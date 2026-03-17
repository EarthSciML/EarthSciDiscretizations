# Operators

## Overview

All finite-volume operators in EarthSciDiscretizations.jl return symbolic `ArrayOp` objects from SymbolicUtils.jl. To obtain numerical results, call `evaluate_arrayop` on the returned `ArrayOp`.

This design allows operators to be composed symbolically before evaluation, enabling future integration with ModelingToolkit.jl equation systems.

## Gradient of a Linear Field

The gradient operator should exactly reproduce the slope of a linear field in computational coordinates.

```@example ops
using EarthSciDiscretizations
using EarthSciDiscretizations: evaluate_arrayop

grid = CubedSphereGrid(8)
Nc = grid.Nc

# Linear field: phi = 2.0 * xi
phi = zeros(6, Nc, Nc)
for p in 1:6, i in 1:Nc, j in 1:Nc
    phi[p, i, j] = 2.0 * grid.ξ_centers[i]
end

ao = fv_gradient_xi(phi, grid)
grad = evaluate_arrayop(ao)
println("Gradient of linear field (slope=2.0):")
println("  Size: $(size(grad))")
println("  Min:  $(minimum(grad))")
println("  Max:  $(maximum(grad))")
println("  All equal to 2.0: $(all(isapprox.(grad, 2.0; rtol=1e-12)))")
```

## Divergence of Zero Flux

The divergence of a zero flux field should be exactly zero.

```@example ops
F_xi = zeros(6, Nc + 1, Nc)
F_eta = zeros(6, Nc, Nc + 1)

ao = fv_divergence(F_xi, F_eta, grid)
div_result = evaluate_arrayop(ao)
println("Divergence of zero flux:")
println("  Size: $(size(div_result))")
println("  Max absolute value: $(maximum(abs.(div_result)))")
```

## Laplacian of a Constant

The Laplacian of a constant field should be zero everywhere.

```@example ops
phi_const = fill(42.0, 6, Nc, Nc)
ao = fv_laplacian(phi_const, grid)
lap = evaluate_arrayop(ao)
println("Laplacian of constant field:")
println("  Size: $(size(lap))")
println("  Max absolute value: $(maximum(abs.(lap)))")
```

## Transport of a Constant Field

Advecting a constant field should produce zero tendency regardless of the velocity field.

```@example ops
q_const = fill(5.0, 6, Nc, Nc)
courant_xi = fill(0.3, 6, Nc, Nc)
courant_eta = fill(-0.2, 6, Nc, Nc)

ao = transport_2d(q_const, courant_xi, courant_eta, grid)
tendency = evaluate_arrayop(ao)
println("Transport tendency of constant field:")
println("  Size: $(size(tendency))")
println("  Max absolute value: $(maximum(abs.(tendency)))")
```

## Operator Summary

| Operator | Input | Output Size | Description |
|:---------|:------|:------------|:------------|
| `fv_divergence` | $F^\xi$ [6, Nc+1, Nc], $F^\eta$ [6, Nc, Nc+1] | [6, Nc, Nc] | Finite-volume divergence |
| `fv_gradient_xi` | $\phi$ [6, Nc, Nc] | [6, Nc-1, Nc] | $\xi$-gradient at interior U-edges |
| `fv_gradient_eta` | $\phi$ [6, Nc, Nc] | [6, Nc, Nc-1] | $\eta$-gradient at interior V-edges |
| `fv_laplacian` | $\phi$ [6, Nc, Nc] | [6, Nc-2, Nc-2] | 5-point Laplacian at interior cells |
| `flux_1d` | $q$, $c$ [6, Nc, Nc] | [6, Nc-2, Nc] or [6, Nc, Nc-2] | 1D Lax-Friedrichs transport |
| `transport_2d` | $q$, $c_\xi$, $c_\eta$ [6, Nc, Nc] | [6, Nc-2, Nc-2] | 2D Lax-Friedrichs transport |
| `ppm_reconstruction` | $q$ [6, Nc, Nc] | Left/Right [6, Nc-4, Nc] | PPM sub-grid reconstruction |
