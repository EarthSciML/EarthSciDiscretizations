# Tutorial: Diffusion on the Sphere

This tutorial simulates diffusion on the cubed sphere using the full
ModelingToolkit integration — you write the PDE, and `discretize` converts
it into an ODE system that can be solved with any SciML solver.

## Problem setup

The diffusion equation on the sphere is:

```math
\frac{\partial u}{\partial t} = \kappa \left(
  \frac{\partial^2 u}{\partial \lambda^2} +
  \frac{\partial^2 u}{\partial \varphi^2}
\right)
```

where ``\lambda`` and ``\varphi`` are longitude and latitude, and ``\kappa``
is the diffusion coefficient. The library automatically maps this notation
to the covariant Laplacian on the cubed-sphere grid.

## Step 1: Define the PDE

```@example tutorial
using EarthSciDiscretizations
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics
using DomainSets
using OrdinaryDiffEqDefault
using SciMLBase

@parameters lon lat
@variables u(..)
Dlon = Differential(lon)
Dlat = Differential(lat)

κ = 0.1  # diffusion coefficient

eq = [D(u(t, lon, lat)) ~ κ * (Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat))))]
bcs = [u(0, lon, lat) ~ exp(-10 * (lon^2 + lat^2))]  # Gaussian blob at (0,0)
domains = [t ∈ Interval(0.0, 0.5),
           lon ∈ Interval(-π, π),
           lat ∈ Interval(-π/2, π/2)]

@named pdesys = PDESystem(eq, bcs, domains, [t, lon, lat], [u(t, lon, lat)])
nothing # hide
```

## Step 2: Discretize and solve

```@example tutorial
Nc = 8  # C8 resolution (6 × 64 = 384 cells)
disc = FVCubedSphere(Nc; R=1.0)

prob = discretize(pdesys, disc)
sol = solve(prob)

println("Retcode: ", sol.retcode)
println("Grid cells: ", 6 * Nc^2)
println("Timesteps: ", length(sol.t))
```

That's the entire workflow: define a PDE with ModelingToolkit, create a
`FVCubedSphere` discretization, and call `discretize`. The library handles
grid construction, spatial derivative replacement, initial condition
projection, and ODE system assembly automatically.

## Step 3: Visualize the result

```@example tutorial
using CairoMakie, GeoMakie

grid = CubedSphereGrid(Nc; R=1.0)

# Extract solution at a given time using symbolic indexing
u_sym = first(@variables u(t)[1:6, 1:Nc, 1:Nc])

function get_snapshot(sol, u_sym, grid, tidx)
    Nc = grid.Nc
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = sol[u_sym[p, i, j]][tidx]
    end
    return q
end

# Interpolate cell-center values to cell corners for vertex coloring
function to_corners(q_panel, Nc)
    c = zeros(Nc + 1, Nc + 1)
    for ii in 1:Nc+1, jj in 1:Nc+1
        n = 0; v = 0.0
        for di in -1:0, dj in -1:0
            ci, cj = ii + di, jj + dj
            if 1 <= ci <= Nc && 1 <= cj <= Nc
                v += q_panel[ci, cj]; n += 1
            end
        end
        c[ii, jj] = v / n
    end
    return c
end

# Build animation of diffusion over time
fig = Figure(size=(900, 500))
ga = GeoAxis(fig[1, 1]; dest="+proj=robin")

q_initial = get_snapshot(sol, u_sym, grid, 1)
color_obs = [Observable(to_corners(q_initial[p, :, :], Nc)) for p in 1:6]

# Use edge coordinates so panels tile seamlessly
for p in 1:6
    lon_corners = zeros(Nc + 1, Nc + 1)
    lat_corners = zeros(Nc + 1, Nc + 1)
    for ii in 1:Nc+1, jj in 1:Nc+1
        lon_corners[ii, jj], lat_corners[ii, jj] =
            gnomonic_to_lonlat(grid.ξ_edges[ii], grid.η_edges[jj], p)
    end
    # Clamp longitudes to [-180, 180] to avoid projection artifacts at map edges
    clamp!(lon_corners, -π, π)
    surface!(ga, rad2deg.(lon_corners), rad2deg.(lat_corners),
             zeros(Nc + 1, Nc + 1); color=color_obs[p], shading=NoShading,
             colormap=:viridis, colorrange=(0, 1))
end
lines!(ga, GeoMakie.coastlines(); color=:black, linewidth=0.5)
Colorbar(fig[1, 2]; colormap=:viridis, colorrange=(0, 1), label="u")

frame_indices = range(1, length(sol.t), length=min(20, length(sol.t))) .|> round .|> Int |> unique

record(fig, "diffusion.gif", frame_indices; framerate=5) do tidx
    q = get_snapshot(sol, u_sym, grid, tidx)
    for p in 1:6
        color_obs[p][] = to_corners(q[p, :, :], Nc)
    end
end
nothing # hide
```

![Diffusion on the sphere](diffusion.gif)

## What happened under the hood

When you call `discretize(pdesys, disc)`:

1. A `CubedSphereGrid` is constructed with the specified resolution
2. For each PDE unknown `u(t, lon, lat)`, a discrete state array `u(t)[1:6, 1:Nc, 1:Nc]` is created
3. The PDE equations are walked symbolically — each spatial `Differential` is replaced with a finite-difference stencil operating on the discrete array
4. Initial conditions from the BCs are evaluated at each grid cell's (lon, lat)
5. The resulting ODE system is compiled by ModelingToolkit and wrapped in an `ODEProblem`

The solver then integrates the ODE system using any algorithm from the
DifferentialEquations.jl ecosystem (Tsit5, Rodas5, etc.).
