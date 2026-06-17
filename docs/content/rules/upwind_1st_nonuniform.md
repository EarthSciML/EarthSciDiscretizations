---
title: "upwind_1st_nonuniform"
slug: "upwind_1st_nonuniform"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx_i)"
applies_to: "grad($u), dim=$x"
rule_path: "discretizations/finite_difference/upwind_1st_nonuniform.json"
description: "Auto-generated catalog entry for the upwind_1st_nonuniform rule (finite_difference, grid family cartesian). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `grad` on the pattern `grad($u), dim=$x` for grid family `cartesian`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$x] = −1/index(dx, $x)·index($u, $x + −1) + 1/index(dx, $x)·index($u, $x)
```

Source: [`discretizations/finite_difference/upwind_1st_nonuniform.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/upwind_1st_nonuniform.json).
