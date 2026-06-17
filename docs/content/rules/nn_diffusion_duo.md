---
title: "nn_diffusion_duo"
slug: "nn_diffusion_duo"
families: "finite_difference"
grid_families: "unstructured"
rule_kinds: "scheme"
accuracy: "O(h^2) on a quasi-uniform triangular icosahedral mesh"
applies_to: "laplacian($u), dim=cell"
rule_path: "discretizations/finite_difference/nn_diffusion_duo.json"
description: "Auto-generated catalog entry for the nn_diffusion_duo rule (finite_difference, grid family unstructured). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `laplacian` on the pattern `laplacian($u), dim=cell` for grid family `unstructured`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out = Σ[k ∈ [0, 2]] ( index(dc_edge, index(edges_on_face, i, k))/(index(dv_edge, index(edges_on_face, i, k))·index(area, i))·(index($u, index(cell_neighbors, i, k)) − index($u, i)) )
```

Source: [`discretizations/finite_difference/nn_diffusion_duo.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/nn_diffusion_duo.json).
