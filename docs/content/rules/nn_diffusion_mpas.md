---
title: "nn_diffusion_mpas"
slug: "nn_diffusion_mpas"
families: "finite_difference"
grid_families: "mpas"
rule_kinds: "scheme"
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

## Discrete operator

For a target cell \(c\) with neighbor list \(\mathcal{N}(c) =
\{c_0, c_1, \dots, c_{n-1}\}\), where \(n = \texttt{n\_edges\_on\_cell}[c]\),
the rule emits

$$\nabla^2 u\big|_c \;\approx\; \sum_{k=0}^{n-1} w_k \bigl(u_{c_k} - u_c\bigr),
\qquad
w_k \;=\; \frac{\ell_v(e_k)}{\ell_c(e_k)\,A_c},$$

where \(e_k = \texttt{edges\_on\_cell}[c, k]\) is the edge separating \(c\) from
neighbor \(c_k\), \(\ell_v(e)\) is its Voronoi (dual) length, \(\ell_c(e)\) the
distance between the two cell centers it bridges, and \(A_c\) the area of cell
\(c\). The rule's two stencil rows materialize the same expression in
two pieces:

| Row | Selector kind | Contributes |
|---|---|---|
| 1 | `reduction` over `cells_on_cell[$target, k]` | \(+w_k\,u_{c_k}\) for each \(k = 0, \dots, n-1\) |
| 2 | `indirect` at `$target` (self) | \(-\bigl(\sum_{k=0}^{n-1} w_k\bigr)\,u_c\) |

Row 2's coefficient is materialized via an `arrayop` reduction so the
diagonal weight is computed from the same `dv_edge / (dc_edge · area_cell)`
expression used in Row 1, guaranteeing the algebraic identity
\(\nabla^2 u(c) = \sum_k w_k\,(u_{c_k} - u_c)\) holds exactly without
relying on a separate analytical simplification.

### Derivation

This stencil is the standard finite-volume Laplacian on an orthogonal
unstructured mesh. Integrating \(\nabla^2 u\) over Voronoi cell \(c\) and
applying the divergence theorem,

$$\int_{c} \nabla^2 u \;dA \;=\; \oint_{\partial c} \nabla u \cdot \mathbf{n}\, ds
\;=\; \sum_{k=0}^{n-1} \int_{e_k} \nabla u \cdot \mathbf{n}\, ds.$$

On an MPAS mesh, the Voronoi edge \(e_k\) is — by construction —
perpendicular to the line connecting cell centers \(c \to c_k\), so the
outward normal flux through \(e_k\) is approximated by the centered
difference

$$\nabla u \cdot \mathbf{n}\big|_{e_k} \;\approx\; \frac{u_{c_k} - u_c}{\ell_c(e_k)},$$

and the integral over the edge contributes

$$\int_{e_k} \nabla u \cdot \mathbf{n}\, ds \;\approx\; \ell_v(e_k)\,\frac{u_{c_k} - u_c}{\ell_c(e_k)}.$$

Summing over edges and dividing by the cell area \(A_c\) (to recover the
pointwise Laplacian from the cell-integrated divergence) yields exactly
the stencil above. This is the diamond-coefficient FV Laplacian; on a
quasi-uniform Voronoi mesh — where edge lengths and cell areas vary
smoothly with \(h\) — the truncation error is \(O(h^2)\), matching the
`accuracy` claim in the rule's frontmatter.

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
| `c_k` | `cells_on_cell` | neighbor cell index, k = 0…n−1 |
| `e_k` | `edges_on_cell` | shared-edge index for each neighbor |
| `n` | `n_edges_on_cell` | per-cell valence (5 or 6 on the icosahedral seed) |

See [the JSON]({{< param repoURL >}}/blob/main/discretizations/finite_difference/nn_diffusion_mpas.json)
for the full coefficient AST and
[`discretizations/SELECTOR_KINDS.md`]({{< param repoURL >}}/blob/main/discretizations/SELECTOR_KINDS.md)
for the `reduction` / `indirect` selector contracts.

## Convergence

<div class="callout callout-pending">
<strong>Pending ESS harness extension.</strong>
The Layer-B convergence harness needs MPAS unstructured support — specifically
the ability to evaluate manufactured solutions on a Voronoi mesh, sample at
cell centers via the MPAS-mesh accessor, and refine across the
<code>x1.642 → x1.2562 → x1.10242</code> sequence. The fixture under
[<code>discretizations/finite_difference/nn_diffusion_mpas/fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_difference/nn_diffusion_mpas/fixtures/convergence)
declares <code>applicable: false</code> until that lands. The structural
claims (rule discoverable, two-row stencil shape, selector kinds, and
coefficient AST containing <code>dv_edge / dc_edge / area_cell</code> plus an
<code>arrayop</code> diagonal-weight reduction) are checked today by
<code>test/test_rule_catalog.jl</code>.
</div>
