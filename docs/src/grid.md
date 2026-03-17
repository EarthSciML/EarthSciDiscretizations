# Cubed-Sphere Grid

## Overview

The cubed-sphere grid projects the six faces of a cube onto the sphere using the gnomonic (central) projection. This avoids the polar singularities of latitude-longitude grids and provides quasi-uniform resolution globally.

Each of the 6 panels is subdivided into $N_c \times N_c$ cells in the local $(\xi, \eta)$ coordinate system, where $\xi, \eta \in [-\pi/4, \pi/4]$.

## Grid Construction


```@example grid
using EarthSciDiscretizations

# Create a C8 grid on the unit sphere
grid = CubedSphereGrid(8; R=1.0)

println("Nc = $(grid.Nc), Ng = $(grid.Ng), R = $(grid.R)")
println("Grid arrays:")
println("  lon/lat: $(size(grid.lon))")
println("  area:    $(size(grid.area))")
println("  dx:      $(size(grid.dx))  (U-edge lengths)")
println("  dy:      $(size(grid.dy))  (V-edge lengths)")
println("  dxi:     $(grid.dξ)")
println("  deta:    $(grid.dη)")
```

## Panel Connectivity

The 6 panels are connected as follows: Panel 1 (front), Panel 2 (right), Panel 3 (top/north pole), Panel 4 (back), Panel 5 (left), Panel 6 (bottom/south pole).


```@example grid
# Show connectivity for Panel 1
for dir in (West, East, South, North)
    nb = PANEL_CONNECTIVITY[1][dir]
    println("Panel 1 $dir -> Panel $(nb.neighbor_panel) ($(nb.neighbor_edge), reverse=$(nb.reverse_index))")
end

# Verify that all connectivity is self-consistent
println("\nConnectivity valid: ", verify_connectivity(8))
```

## Cell Area Statistics

```@example grid
using Statistics

panel_areas = [sum(grid.area[p, :, :]) for p in 1:6]
cell_areas = vec(grid.area)

println("Total area: $(sum(grid.area))")
println("Expected (4pi): $(4pi)")
println("Panel areas: ", round.(panel_areas; digits=6))
println("Min cell area: $(minimum(cell_areas))")
println("Max cell area: $(maximum(cell_areas))")
println("Max/Min ratio: $(maximum(cell_areas) / minimum(cell_areas))")
println("Std/Mean: $(std(cell_areas) / mean(cell_areas))")
```

## Projection Functions

