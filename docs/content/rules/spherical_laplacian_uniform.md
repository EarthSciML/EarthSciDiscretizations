---
title: "spherical_laplacian_uniform"
slug: "spherical_laplacian_uniform"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dr^2)"
applies_to: "spherical_laplacian($u), dim=$r"
rule_path: "discretizations/finite_difference/spherical_laplacian_uniform.json"
description: "Auto-generated catalog entry for the spherical_laplacian_uniform rule (finite_difference, grid family cartesian). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `spherical_laplacian` on the pattern `spherical_laplacian($u), dim=$r` for grid family `cartesian`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$r] = ((0.5·(index(r, $r) + index(r, $r + 1)))²·(index($u, $r + 1) − index($u, $r)) − (0.5·(index(r, $r − 1) + index(r, $r)))²·(index($u, $r) − index($u, $r − 1)))/(index(r, $r)²·dr·dr)
```

Source: [`discretizations/finite_difference/spherical_laplacian_uniform.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/spherical_laplacian_uniform.json).
