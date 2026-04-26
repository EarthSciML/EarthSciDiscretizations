---
title: "Rule × grid matrix"
description: "Auto-generated coverage matrix: which rules dispatch on which grid families, and the convergence-order observed by each Layer-B fixture."
---

The matrix below is regenerated from
[`discretizations/`]({{< param repoURL >}}/tree/main/discretizations) and the
fixture trees under
[`tests/conformance/`]({{< param repoURL >}}/tree/main/tests/conformance) by
[`tools/render_rule_matrix.py`]({{< param repoURL >}}/tree/main/tools/render_rule_matrix.py).
**Drop a new rule or a new grid into the catalog and CI updates this view —
no hardcoded cells.** Use the matrix to answer two questions:

1. *(user)* "Does rule X work on grid Y?" — applicable cells are coloured;
   the number is the floor on observed convergence order baked into the
   rule's Layer-B fixture (`expected_min_order` in
   `discretizations/<family>/<rule>/fixtures/convergence/expected.esm`).
2. *(contributor)* "Where are the holes?" — `missing` cells flag rule × grid
   pairs the rule declares it serves but where the canonical or convergence
   fixture has not landed yet.

Cells legend:

- **`order` (green)** — applicable, convergence fixture is running, value is
  the fixture's `expected_min_order` floor.
- **`canonical` (blue)** — applicable, cross-binding canonical fixture
  exists but no convergence fixture yet.
- **`missing` (amber)** — applicable, no fixture yet (or fixture is
  structurally skipped pending an ESS harness extension). Visible dispatch
  hole.
- **`n/a` (muted)** — rule does not declare this grid family in its
  `applies_to` block.
