# EarthSciDiscretizations.jl

**EarthSciDiscretizations.jl** provides finite-volume (FV) discretization of partial differential equations on the cubed-sphere grid. It is designed for global atmospheric and oceanic modeling, where the cubed-sphere avoids the polar singularities of latitude-longitude grids while providing quasi-uniform resolution.

The transport algorithms are based on the Lin-Rood (1996) dimensionally-split scheme with Piecewise Parabolic Method (PPM) reconstruction (Colella & Woodward 1984), adapted to the cubed-sphere geometry following Putman & Lin (2007). Vertical remapping follows Lin (2004). See [Finite-Volume Method](@ref) for full references.

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
- [Tutorial: Authoring a rule](@ref) -- end-to-end walkthrough for writing a new discretization rule in the closed-AST lowering pattern
