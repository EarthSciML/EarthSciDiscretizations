# EarthSciDiscretizations.jl

**EarthSciDiscretizations.jl** provides finite-volume (FV) discretization of partial differential equations on the cubed-sphere grid. It is designed for global atmospheric and oceanic modeling, where the cubed-sphere avoids the polar singularities of latitude-longitude grids while providing quasi-uniform resolution.

## Features

- **Cubed-sphere grid construction** with gnomonic equidistant projection, automatic metric tensor computation, and panel connectivity
- **C-grid staggering** for scalar fields (cell centers), velocity components (U/V edges), and corner quantities
- **Ghost cell management** for seamless inter-panel communication
- **Finite-volume operators** implemented as symbolic `ArrayOp` expressions: divergence, gradient, Laplacian, PPM reconstruction, 1D flux, and 2D transport
- **Discretization pipeline** with `FVCubedSphere` specification, operator registry, and initial condition projection

## Quick Start

```@example quickstart
using EarthSciDiscretizations

# Create a C16 cubed-sphere grid (16x16 cells per panel, unit sphere)
grid = CubedSphereGrid(16; R=1.0)

println("Grid: $(grid.Nc)x$(grid.Nc) cells per panel, 6 panels")
println("Total cells: $(6 * grid.Nc^2)")
println("Total area: $(total_area(grid))")
println("Expected (4pi): $(4pi)")
println("Relative error: $(abs(total_area(grid) - 4pi) / 4pi)")
```

## Package Contents

- [Finite-Volume Method](@ref) -- mathematical foundations and operator formulas
- [Cubed-Sphere Grid](@ref) -- grid construction and geometry
- [Operators](@ref) -- using FV operators on data
- [Tutorial: Diffusion on the Sphere](@ref) -- full diffusion simulation example
