---
title: "column_nonuniform_1d_tracer"
description: "grids/column_nonuniform_1d with a trailing TRACER axis: fields shaped [lev, tracer] are NTRACER transported constituents (chemical species, aerosol sections) that share ONE column's geometry and coefficients."
esd_kind: "grid"
families: ["column", "finite_volume"]
ops: ["interface_flux", "varcoeff_laplacian"]
bcs: ["prescribed_flux", "supplied_face_flux"]
tags: ["esd:grid", "family:column", "ndim:2", "spacing:nonuniform", "axis:lev", "axis:tracer", "esd:stencil", "family:finite_volume", "grid:column_nonuniform_1d_tracer", "op:interface_flux", "order:2", "stencil_width:2", "op:varcoeff_laplacian", "stencil_width:3", "esd:rule", "axes:lev,tracer", "bc:prescribed_flux", "spec:9.6.8", "bc:supplied_face_flux"]
source: "grids/column_nonuniform_1d_tracer/grid.typ"
---

grids/column_nonuniform_1d with a trailing TRACER axis: fields shaped [lev, tracer] are NTRACER transported constituents (chemical species, aerosol sections) that share ONE column's geometry and coefficients. The vertical axes are exactly the column grid's -- `lev` (NLEV layers, layer 1 the surface) and `lev_nodes` (NLEV+1 interfaces) -- declared here -- with the NLEV metaparameter -- verbatim as on that grid, so they deep-equal and dedupe against grids/column_nonuniform_1d when a model imports from both grids; `tracer` (NTRACER members) is the only new axis. The geometry is the column grid's, unchanged: a consumer supplies the interface array `ze` [lev_nodes] and derives `zc` and `dz` [lev] by applying column_nonuniform_1d's `column_cell_center` / `column_cell_thickness` (import that grid file, or its mesh fixture, beside this grid's rules); this grid adds no geometry, so the per-tracer operators read the same `zc`/`dz` by bare name. Its one template, `tracer_column_total`, is the per-tracer column total sum_k u[k,s] dz[k] that the flux-form operators here conserve tracer by tracer. It is a separate grid rather than an extra index set on column_nonuniform_1d because the grid-first layout (AGENTS.md section 3) puts an operator that takes a different form on a different field shape under its own grid, and rule match-scoping (esm-spec 9.6.1 `where` shapes) already tells a [lev] field from a [lev, tracer] one by axis names, so the two grids' rules are disjoint and a model may import both. A consumer whose species axis has another name imports this grid's rules with `rename: {tracer: <its axis>}` (esm-spec 9.7.7) and binds NTRACER to that axis's size; NLEV is bound as on the column grid. Motivation: one mount of a column tracer-mixing component for every species of an atmospheric-chemistry mechanism (WRF-Chem's vertmx is called once per species with the same exchange coefficients).

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NLEV` | integer | `64` | Number of layers (cells) along lev; the interface axis lev_nodes has NLEV+1 members. Minimum 2 (see the grid description). |
| `NTRACER` | integer | `1` | Number of tracers on the trailing axis. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Index sets

| Name | Kind | Size |
|---|---|---|
| `lev` | interval | \(\mathit{NLEV}\) |
| `lev_nodes` | interval | \(\mathit{NLEV} + 1\) |
| `tracer` | interval | \(\mathit{NTRACER}\) |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

