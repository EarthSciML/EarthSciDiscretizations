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


### Gradient

The discrete gradient is computed at edge locations. The $\xi$-component at interior U-edges:

```math
\left(\frac{\partial \phi}{\partial \xi}\right)_{p,i,j} = \frac{\phi_{p,i+1,j} - \phi_{p,i,j}}{\Delta \xi}
```

The $\eta$-component at interior V-edges:

```math
\left(\frac{\partial \phi}{\partial \eta}\right)_{p,i,j} = \frac{\phi_{p,i,j+1} - \phi_{p,i,j}}{\Delta \eta}
```


### Laplacian

The full covariant Laplacian on the cubed sphere is:

```math
\nabla^2 \phi = \frac{1}{J} \left[ \frac{\partial}{\partial \xi}\left(J\, g^{\xi\xi} \frac{\partial \phi}{\partial \xi} + J\, g^{\xi\eta} \frac{\partial \phi}{\partial \eta}\right) + \frac{\partial}{\partial \eta}\left(J\, g^{\xi\eta} \frac{\partial \phi}{\partial \xi} + J\, g^{\eta\eta} \frac{\partial \phi}{\partial \eta}\right) \right]
```

where $J$ is the Jacobian determinant and $g^{\alpha\beta}$ are the inverse metric tensor components.

The discrete operator uses a 5-point orthogonal stencil plus a cross-metric correction:

```math
\nabla^2 \phi_{p,i,j} \approx \frac{1}{A_{p,i,j}} \left[ \text{(face gradient} \times \text{edge length)} \right] + 2\, g^{\xi\eta} \frac{\partial^2 \phi}{\partial \xi \, \partial \eta} + \frac{1}{J}\left[\frac{\partial(J\, g^{\xi\eta})}{\partial \xi} \frac{\partial \phi}{\partial \eta} + \frac{\partial(J\, g^{\xi\eta})}{\partial \eta} \frac{\partial \phi}{\partial \xi}\right]
```

The orthogonal part uses physical center-to-center distances and edge lengths. The cross-metric $g^{\xi\eta}$ correction uses a 4-point cross-stencil for the mixed second derivative and includes first-derivative correction terms arising from the spatial variation of $J \cdot g^{\xi\eta}$ to account for the non-orthogonality of the gnomonic cubed-sphere grid.


### Transport

The 1D and 2D transport operators use a Lax-Friedrichs (upwind) flux splitting. The numerical flux at edge $(i+1/2)$ is:

```math
F_{i+1/2} = \frac{c_{i+1/2} (q_i + q_{i+1})}{2} - \frac{|c_{i+1/2}| (q_{i+1} - q_i)}{2}
```

where $c$ is the Courant number and $q$ is the transported scalar. The tendency is:

```math
\frac{\partial q}{\partial t} = -\frac{F_{i+1/2} - F_{i-1/2}}{\Delta \xi}
```


### PPM Reconstruction

The Piecewise Parabolic Method (PPM) provides higher-order sub-grid reconstruction for transport. It computes left and right interface values using a 4th-order interpolation formula with monotonicity limiting (Colella & Woodward, 1984).


### Vertical Remapping

Vertical remapping uses PPM with Colella-Woodward (1984) monotonicity limiting to conservatively remap quantities between different vertical layer structures. The remapping preserves column-integrated mass ($\sum q \cdot \Delta p$) exactly. See `vertical_remap`.


## C-Grid Staggering

The Arakawa C-grid staggering places different variables at different locations within each cell:

| Location | Symbol | Grid Size | Description |
|:---------|:-------|:----------|:------------|
| `CellCenter` | $(i, j)$ | $(N_c, N_c)$ | Scalar fields (tracer, pressure, temperature) |
| `UEdge` | $(i+1/2, j)$ | $(N_c+1, N_c)$ | Normal velocity component in $\xi$-direction |
| `VEdge` | $(i, j+1/2)$ | $(N_c, N_c+1)$ | Normal velocity component in $\eta$-direction |
| `Corner` | $(i+1/2, j+1/2)$ | $(N_c+1, N_c+1)$ | Vorticity, stream function |


## Ghost Cells

Inter-panel communication is handled through ghost cells. Each panel is padded with $N_g$ ghost layers on each side, filled from neighboring panels using the connectivity table and index transformations.


## ArrayOp Utilities

All operators return symbolic `ArrayOp` objects. Use `evaluate_arrayop` to obtain numerical results.


## References

The finite-volume methods implemented in this package are based on the following foundational works:

- Lin, S.-J. and R. B. Rood (1996). "Multidimensional Flux-Form Semi-Lagrangian Transport Schemes." *Monthly Weather Review*, 124(9), 2046--2070. — Dimensionally-split transport algorithm.
- Colella, P. and P. R. Woodward (1984). "The Piecewise Parabolic Method (PPM) for gas-dynamical simulations." *Journal of Computational Physics*, 54(1), 174--201. — PPM reconstruction and monotonicity limiter.
- Lin, S.-J. (2004). "A 'Vertically Lagrangian' Finite-Volume Dynamical Core for Global Models." *Monthly Weather Review*, 132(10), 2293--2307. — Vertically Lagrangian FV framework and conservative vertical remapping.
- Putman, W. M. and S.-J. Lin (2007). "Finite-volume transport on various cubed-sphere grids." *Journal of Computational Physics*, 227(1), 55--78. — FV transport on cubed-sphere grids.
- Ronchi, C., R. Iacono, and P. S. Paolucci (1996). "The 'Cubed Sphere': A New Method for the Solution of Partial Differential Equations in Spherical Geometry." *Journal of Computational Physics*, 124(1), 93--114. — Gnomonic cubed-sphere projection and metric tensors.

