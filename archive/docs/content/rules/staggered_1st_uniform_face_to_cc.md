---
title: "staggered_1st_uniform_face_to_cc"
slug: "staggered_1st_uniform_face_to_cc"
families: "finite_difference"
grid_families: "arakawa"
rule_kinds: "scheme"
accuracy: "O(h^2) (two-point centered difference at the half-offset cell center)"
applies_to: "grad($u)"
rule_path: "discretizations/finite_difference/staggered_1st_uniform.json"
description: "Auto-generated catalog entry for the staggered_1st_uniform_face_to_cc rule (finite_difference, grid family arakawa). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `grad` on the pattern `grad($u)` for grid family `arakawa`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$x] = −1/h·index($u, $x) + 1/h·index($u, $x + 1)
```

Source: [`discretizations/finite_difference/staggered_1st_uniform.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/staggered_1st_uniform.json).
