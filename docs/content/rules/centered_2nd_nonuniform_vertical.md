---
title: "centered_2nd_nonuniform_vertical"
slug: "centered_2nd_nonuniform_vertical"
families: "finite_difference"
grid_families: "vertical"
rule_kinds: "scheme"
accuracy: "O(dz_k^2)"
applies_to: "grad($u), dim=$k"
rule_path: "discretizations/finite_difference/centered_2nd_nonuniform_vertical.json"
description: "Auto-generated catalog entry for the centered_2nd_nonuniform_vertical rule (finite_difference, grid family vertical). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `grad` on the pattern `grad($u), dim=$k` for grid family `vertical`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
(index($u, $k + 1) − index($u, $k − 1))/(2·index(dz, $k))
```

Source: [`discretizations/finite_difference/centered_2nd_nonuniform_vertical.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/centered_2nd_nonuniform_vertical.json).
