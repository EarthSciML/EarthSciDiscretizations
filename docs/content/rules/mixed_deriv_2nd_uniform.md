---
title: "mixed_deriv_2nd_uniform"
slug: "mixed_deriv_2nd_uniform"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx^2)"
applies_to: "grad({'op': 'grad', 'args': ['$u'], 'dim': '$y'}), dim=$x"
rule_path: "discretizations/finite_difference/mixed_deriv_2nd_uniform.json"
description: "Auto-generated catalog entry for the mixed_deriv_2nd_uniform rule (finite_difference, grid family cartesian). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `grad` on the pattern `grad({'op': 'grad', 'args': ['$u'], 'dim': '$y'}), dim=$x` for grid family `cartesian`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$x, $y] = (index($u, $x + 1, $y + 1) − index($u, $x + 1, $y − 1) − (index($u, $x − 1, $y + 1) − index($u, $x − 1, $y − 1)))/(4·dx·dy)
```

Source: [`discretizations/finite_difference/mixed_deriv_2nd_uniform.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/mixed_deriv_2nd_uniform.json).
