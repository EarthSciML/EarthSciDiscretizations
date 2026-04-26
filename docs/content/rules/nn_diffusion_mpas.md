---
title: "nn_diffusion_mpas"
slug: "nn_diffusion_mpas"
family: "finite_difference"
grid_family: "mpas"
rule_kind: "scheme"
accuracy: "O(h²) on quasi-uniform Voronoi"
applies_to: "laplacian(u), dim=cell"
rule_path: "discretizations/finite_difference/nn_diffusion_mpas.json"
description: "Nearest-neighbor diffusion on an MPAS Voronoi mesh, summed over edges_on_cell with dv_edge / (dc_edge · area_cell) coefficients."
tags: ["finite-difference", "mpas", "voronoi", "diffusion", "unstructured"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/nn_diffusion_mpas-stencil.png"
       alt="Voronoi neighbor stencil — central cell with 6 hexagonal neighbors">
  <figcaption>Illustrative stencil over a hexagonal Voronoi cell. The actual
  number of neighbors per cell is <code>n_edges_on_cell[c]</code> — 5 for the
  twelve pentagons inherited from the icosahedral seed, 6 elsewhere.</figcaption>
</figure>

## Coefficients

The Laplacian is approximated by summing flux contributions across each edge
of the target cell:

$$\nabla^2 u\big|_c \approx \sum_{k=0}^{n_{edges}(c)-1} \frac{\ell_v(e_k)}{\ell_c(e_k)\,A_c} \bigl(u_{c_k} - u_c\bigr),$$

with the edge metric quantities loaded from the MPAS mesh tables:

| Symbol | MPAS table | Meaning |
|---|---|---|
| `ℓ_v` | `dv_edge` | Voronoi edge length (between dual vertices) |
| `ℓ_c` | `dc_edge` | dual-edge length (between cell centers) |
| `A_c` | `area_cell` | cell area |

The rule uses two stencil entries: a `reduction` selector that sweeps the
neighbor index `k` over `cells_on_cell[$target, k]` and emits the positive
neighbor contribution, and an `indirect` selector that emits the negated
sum at the central cell. See `discretizations/SELECTOR_KINDS.md`.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B convergence harness needs MPAS unstructured support — specifically
the ability to evaluate manufactured solutions on a Voronoi mesh, sample at
cell centers via the MPAS-mesh accessor, and refine across the
<code>x1.642 → x1.2562 → x1.10242</code> sequence. The fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/nn_diffusion_mpas/fixtures/convergence)
declares <code>applicable: false</code> until that lands.
</div>
