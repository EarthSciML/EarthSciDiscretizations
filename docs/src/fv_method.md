# Finite-Volume Method

## Overview

EarthSciDiscretizations.jl uses the finite-volume (FV) method to discretize PDEs on the cubed-sphere grid. The FV method is based on the integral form of conservation laws, making it naturally conservative and well-suited for geophysical fluid dynamics.

## The Divergence Theorem

The foundation of the FV method is the divergence theorem. For a vector field $\mathbf{F}$ over a cell with area $A$ and boundary $\partial A$:

```math
\int_A \nabla \cdot \mathbf{F} \, dA = \oint_{\partial A} \mathbf{F} \cdot \hat{n} \, ds
```

In the discrete setting, each cell is a quadrilateral on the cubed-sphere surface. The integral of the divergence is approximated by the sum of fluxes through the cell edges.

## Discrete Operators

### Divergence

The discrete divergence at cell center $(i, j)$ on panel $p$ is:

```math
(\nabla \cdot \mathbf{F})_{p,i,j} = \frac{1}{A_{p,i,j}} \left[ F^\xi_{p,i+1,j} \, \Delta x_{p,i+1,j} - F^\xi_{p,i,j} \, \Delta x_{p,i,j} + F^\eta_{p,i,j+1} \, \Delta y_{p,i,j+1} - F^\eta_{p,i,j} \, \Delta y_{p,i,j} \right]
```

where $F^\xi$ and $F^\eta$ are the contravariant flux components at U-edges and V-edges respectively, $\Delta x$ and $\Delta y$ are edge lengths, and $A$ is the cell area.

```@docs
fv_divergence
```

### Gradient

The discrete gradient is computed at edge locations. The $\xi$-component at interior U-edges:

```math
\left(\frac{\partial \phi}{\partial \xi}\right)_{p,i,j} = \frac{\phi_{p,i+1,j} - \phi_{p,i,j}}{\Delta \xi}
```

The $\eta$-component at interior V-edges:

```math
\left(\frac{\partial \phi}{\partial \eta}\right)_{p,i,j} = \frac{\phi_{p,i,j+1} - \phi_{p,i,j}}{\Delta \eta}
```

```@docs
fv_gradient_xi
fv_gradient_eta
```

### Laplacian

The discrete Laplacian uses a 5-point stencil with metric corrections:

```math
\nabla^2 \phi_{p,i,j} = \frac{1}{A_{p,i,j}} \left[ \frac{\phi_{p,i+1,j} - \phi_{p,i,j}}{\Delta \xi} \Delta x_{p,i+1,j} - \frac{\phi_{p,i,j} - \phi_{p,i-1,j}}{\Delta \xi} \Delta x_{p,i,j} + \frac{\phi_{p,i,j+1} - \phi_{p,i,j}}{\Delta \eta} \Delta y_{p,i,j+1} - \frac{\phi_{p,i,j} - \phi_{p,i,j-1}}{\Delta \eta} \Delta y_{p,i,j} \right]
```

The edge lengths $\Delta x$ and $\Delta y$ account for the non-uniform metric of the gnomonic projection.

```@docs
fv_laplacian
```

### Transport

The 1D and 2D transport operators use a Lax-Friedrichs (upwind) flux splitting. The numerical flux at edge $(i+1/2)$ is:

```math
F_{i+1/2} = \frac{c_{i+1/2} (q_i + q_{i+1})}{2} - \frac{|c_{i+1/2}| (q_{i+1} - q_i)}{2}
```

where $c$ is the Courant number and $q$ is the transported scalar. The tendency is:

```math
\frac{\partial q}{\partial t} = -\frac{F_{i+1/2} - F_{i-1/2}}{\Delta \xi}
```

```@docs
flux_1d
transport_2d
```

### PPM Reconstruction

The Piecewise Parabolic Method (PPM) provides higher-order sub-grid reconstruction for transport. It computes left and right interface values using a 4th-order interpolation formula with monotonicity limiting (Colella & Woodward, 1984).

```@docs
ppm_reconstruction
ppm_reconstruction!
```

### Vertical Remapping

Vertical remapping using PPM is planned but not yet implemented.

```@docs
vertical_remap_tendency
```

## C-Grid Staggering

The Arakawa C-grid staggering places different variables at different locations within each cell:

| Location | Symbol | Grid Size | Description |
|:---------|:-------|:----------|:------------|
| `CellCenter` | $(i, j)$ | $(N_c, N_c)$ | Scalar fields (tracer, pressure, temperature) |
| `UEdge` | $(i+1/2, j)$ | $(N_c+1, N_c)$ | Normal velocity component in $\xi$-direction |
| `VEdge` | $(i, j+1/2)$ | $(N_c, N_c+1)$ | Normal velocity component in $\eta$-direction |
| `Corner` | $(i+1/2, j+1/2)$ | $(N_c+1, N_c+1)$ | Vorticity, stream function |

```@docs
VarLocation
grid_size
full_array_size
ghost_array_size
```

## Ghost Cells

Inter-panel communication is handled through ghost cells. Each panel is padded with $N_g$ ghost layers on each side, filled from neighboring panels using the connectivity table and index transformations.

```@docs
fill_ghost_cells!
extend_with_ghosts
```

## ArrayOp Utilities

All operators return symbolic `ArrayOp` objects. Use `evaluate_arrayop` to obtain numerical results.

```@docs
const_wrap
get_idx_vars
make_arrayop
evaluate_arrayop
```
