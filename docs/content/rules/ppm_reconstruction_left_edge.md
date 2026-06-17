---
title: "ppm_reconstruction_left_edge"
slug: "ppm_reconstruction_left_edge"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx^4)"
applies_to: "q_left_edge($q), dim=$x"
rule_path: "discretizations/finite_volume/ppm_reconstruction_left_edge.json"
description: "Auto-generated catalog entry for the ppm_reconstruction_left_edge rule (finite_volume, grid family cartesian). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `q_left_edge` on the pattern `q_left_edge($q), dim=$x` for grid family `cartesian`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$x] = (−1·index($q, $x − 2) + 7·index($q, $x − 1))/12
```

Source: [`discretizations/finite_volume/ppm_reconstruction_left_edge.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_reconstruction_left_edge.json).
