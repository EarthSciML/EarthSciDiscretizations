---
title: "EarthSciDiscretizations"
description: "Catalog of declarative discretization rules and grid families used by ESD."
---

# EarthSciDiscretizations catalog

This site is auto-generated from the rule files under
[`discretizations/`]({{< param repoURL >}}/tree/main/discretizations) and the
grid family sources under
[`src/grids/`]({{< param repoURL >}}/tree/main/src/grids).

Each **grid family** page describes the topology of a supported grid
(connectivity, coordinates, metric tensor) and shows a visualization of a
typical configuration. Each **rule** page renders the stencil pattern and
coefficient diagram, plus a convergence plot for rules whose Layer-B fixtures
are currently producing.

- **[Grid families →]({{< ref "/grids" >}})**
- **[Discretization rules →]({{< ref "/rules" >}})**
- **[Tutorials →]({{< ref "/tutorials" >}})** — contributor walkthroughs,
  starting with [adding a new rule]({{< ref "/tutorials/add-a-rule" >}}).
- **[Rule × grid matrix →]({{< ref "/matrix" >}})** — what dispatches where,
  with observed convergence orders and visible holes in fixture coverage.
- Browse [by family]({{< ref "/families" >}}) (`finite_difference`,
  `finite_volume`, …)
- Browse [by grid family]({{< ref "/grid_families" >}})

The source of truth is the rule file (JSON) and the grid trait
implementation (Julia). Drop a new rule under `discretizations/<family>/`
and CI regenerates this catalog.
