---
title: "weno5_grad"
slug: "weno5_grad"
families: "finite_difference"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(dx^5)"
applies_to: "grad($u), dim=$x"
rule_path: "discretizations/finite_difference/weno5_grad.json"
description: "Auto-generated catalog entry for the weno5_grad rule (finite_difference, grid family cartesian). The discrete operator is pretty-printed directly from the rule's replacement AST."
---

<div class="callout">
This page is <strong>auto-generated</strong> from the rule's JSON by <code>tools/render_doc_plots.py</code>. The discrete operator below is pretty-printed directly from the rule's replacement AST using EarthSciSerialization (<code>earthsci_toolkit</code>). Rules with a hand-authored page carry richer prose plus stencil and convergence diagrams; see the curated entries linked from the <a href="{{< ref "/rules" >}}">rule index</a>.
</div>

## Pattern

Matches the §4.2 operator `grad` on the pattern `grad($u), dim=$x` for grid family `cartesian`.

## Discrete operator

Pretty-printed from the rule's `replacement` AST (`$`-prefixed names are pattern variables bound by `applies_to`; grid metrics such as `dx`/`dz` are resolved by the grid):

```text
out[$x] = ((1/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 2) + −2·index($u, $x − 1))² + 1/4·(index($u, $x − 2) + −4·index($u, $x − 1))²)²·(11/6·index($u, $x) + −7/6·index($u, $x − 1)) + 6/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 1) + −2·index($u, $x))² + 1/4·(index($u, $x − 1) + −1·index($u, $x + 1))²)²·(−1/6·index($u, $x − 1) + 5/6·index($u, $x)))/(1/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 2) + −2·index($u, $x − 1))² + 1/4·(index($u, $x − 2) + −4·index($u, $x − 1))²)² + 6/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 1) + −2·index($u, $x))² + 1/4·(index($u, $x − 1) + −1·index($u, $x + 1))²)²) − (1/10/(1.0×10⁻⁶ + 13/12·(index($u, $x) + −2·index($u, $x + 1))² + 1/4·(3·index($u, $x) + −4·index($u, $x + 1))²)²·(11/6·index($u, $x) + −7/6·index($u, $x + 1)) + 6/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 1) + −2·index($u, $x))² + 1/4·(index($u, $x − 1) + −1·index($u, $x + 1))²)²·(1/3·index($u, $x − 1) + 5/6·index($u, $x)))/(1/10/(1.0×10⁻⁶ + 13/12·(index($u, $x) + −2·index($u, $x + 1))² + 1/4·(3·index($u, $x) + −4·index($u, $x + 1))²)² + 6/10/(1.0×10⁻⁶ + 13/12·(index($u, $x − 1) + −2·index($u, $x))² + 1/4·(index($u, $x − 1) + −1·index($u, $x + 1))²)²))/dx
```

Source: [`discretizations/finite_difference/weno5_grad.json`]({{< param repoURL >}}/blob/main/discretizations/finite_difference/weno5_grad.json).
