---
title: "lax_friedrichs_flux_cubed_sphere"
slug: "lax_friedrichs_flux_cubed_sphere"
families: "finite_volume"
grid_families: "cubed_sphere"
rule_kinds: "scheme"
accuracy: "O(h)"
applies_to: "flux(q, c, dim ∈ {xi, eta})"
rule_path: "discretizations/finite_volume/lax_friedrichs_flux_cubed_sphere_xi.json"
description: "Lax-Friedrichs face flux on the gnomonic cubed sphere — sibling-pair (ξ at u_edge, η at v_edge) wrappers for the Cartesian LF flux algebra with face-staggered Courant bindings."
tags: ["finite-volume", "cubed-sphere", "lax-friedrichs", "transport", "face-stagger"]
---

This page covers the sibling-file pair
[`lax_friedrichs_flux_cubed_sphere_xi.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux_cubed_sphere_xi.json)
and
[`lax_friedrichs_flux_cubed_sphere_eta.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux_cubed_sphere_eta.json).
The two rules are coordinate-mirror images of each other under
`(ξ, west, east) ↔ (η, south, north)`; the discussion here uses the
ξ-direction rule unless otherwise noted.

## Continuous form

For the linear-advection flux `f(q) = c·q` on a per-face Courant `c`, the
Lax & Friedrichs (1954) numerical flux at an interface `i+1/2` is

$$ F_{i+1/2} \;=\; \tfrac{c}{2}\,(q_i + q_{i+1}) \;-\; \tfrac{|c|}{2}\,(q_{i+1} - q_i). $$

Algebraically this reduces to first-order upwinding,

$$ F_{i+1/2} \;=\; \tfrac{c + |c|}{2}\,q_i \;+\; \tfrac{c - |c|}{2}\,q_{i+1} \;=\; \max(c, 0)\,q_i + \min(c, 0)\,q_{i+1}, $$

with no caller-side branching on `sign(c)` — the upwind cell selection is
encoded directly through the `abs` op. The Cartesian core of this identity
lives in
[`lax_friedrichs_flux.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux.json);
the cubed-sphere wrappers swap selector kind / stagger but keep the
algebra coefficient-for-coefficient.

## Discrete operator

The ξ-direction rule emits at `u_edge` stagger and reads two cell-center
samples flanking the face (offsets `{-1, 0}` along ξ at `cell_center`)
plus the face's own Courant (offset `0` along ξ at `u_edge`):

| offset along `ξ` | stagger | coefficient | role |
|---|---|---|---|
| `-1` | `cell_center` | `(c + abs(c)) / 2` | west cell `q_west`; contributes when `c > 0` |
| ` 0` | `cell_center` | `(c - abs(c)) / 2` | east cell `q_east`; contributes when `c < 0` |
| ` 0` | `u_edge`      | (binding `$c`)    | face-staggered Courant precomputed by the cubed_sphere accessor |

The η sibling has the analogous `v_edge` form with offsets along η:

| offset along `η` | stagger | coefficient | role |
|---|---|---|---|
| `-1` | `cell_center` | `(c + abs(c)) / 2` | south cell `q_south` |
| ` 0` | `cell_center` | `(c - abs(c)) / 2` | north cell `q_north` |
| ` 0` | `v_edge`      | (binding `$c`)    | face Courant |

## Cubed-sphere face-stagger schema

These rules are the first cubed-sphere flux schemes to combine
[SELECTOR_KINDS.md decision #17]({{< param repoURL >}}/blob/main/discretizations/SELECTOR_KINDS.md)
(per-axis stagger field on the `cubed_sphere` selector — `cell_center`,
`u_edge`, `v_edge`, `corner`) with the existing per-cell metric-binding
contract from
[`covariant_laplacian_cubed_sphere`]({{< relref "covariant_laplacian_cubed_sphere" >}})
and
[`ppm_edge_cubed_sphere`]({{< relref "ppm_edge_cubed_sphere" >}}).

- Selectors **do not** carry a `panel` field (decision #13). Cross-panel
  ghost extension at the western (i=1) and eastern (i=Nc+1) panel
  boundaries — and the analogous south/north η boundaries — is resolved
  by the cubed_sphere grid accessor
  ([`src/grids/panel_connectivity.jl`]({{< param repoURL >}}/blob/main/src/grids/panel_connectivity.jl))
  via `extend_with_ghosts` for the scalar `$q`.
- The face-staggered Courant `$c` is consumed via a `reads` block at
  `u_edge` / `v_edge` stagger — the same pattern used by
  [`fv3_sinsg_flux_xi`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/fv3_sinsg_flux_xi.json)
  for its face velocity. Panel-boundary distance handling
  (`grid.dist_xi_bnd` / `grid.dist_eta_bnd` at the seam, vs `grid.dist_xi`
  / `grid.dist_eta` in the interior) is encapsulated in the imperative
  Courant precomputation `_get_courant_xi` / `_get_courant_eta` /
  `compute_courant_numbers_arrayop` in
  [`src/operators/flux_1d.jl`]({{< param repoURL >}}/blob/main/src/operators/flux_1d.jl);
  by the time the rule sees `$c` it is a single face-staggered array
  with boundary distances already factored in.

A `boundary_policy` schema feature that lets the rule reference
`grid.dist_xi_bnd` / `grid.dist_eta_bnd` directly is a tracked follow-up
off this bead (dsc-0fd) — see the rule's `schema_notes`.

## Composition

Caller workflow (mirrors `transport_2d_linrood!` / `flux_1d` in the
imperative reference):

1. Precompute face-staggered Courants
   `c[p, i, j] = u_face[p, i, j] · dt / dist[p, i, j]` via
   `compute_courant_numbers_arrayop` ([`src/operators/flux_1d.jl`]({{< param repoURL >}}/blob/main/src/operators/flux_1d.jl)).
2. Apply this rule at every ξ-face on every panel to produce
   `F_xi[p, i, j]` at `u_edge` stagger; apply the η sibling for
   `F_eta[p, i, j]` at `v_edge`.
3. Match same-direction panel-boundary fluxes via
   `_match_boundary_fluxes_xi!` / `_match_boundary_fluxes_eta!` (mass
   averaging across the seam). For rotated ξ↔η panel connections, run
   `_match_rotated_boundary_fluxes!` after both `F_xi` and `F_eta` are
   computed. These conservation post-processing steps are not yet
   expressible in the rule schema (follow-up shared with
   [`ppm_edge_cubed_sphere`]({{< relref "ppm_edge_cubed_sphere" >}})'s
   `monotonicity_constraint`).
4. Convert face fluxes to a cell-centered tendency via the FV divergence:

$$ T_{p,i,j} \;=\; -\frac{F^\xi_{p,i+1,j}\,dx_{p,i+1,j} - F^\xi_{p,i,j}\,dx_{p,i,j}}{\mathrm{area}_{p,i,j}} \;-\; \frac{F^\eta_{p,i,j+1}\,dy_{p,i,j+1} - F^\eta_{p,i,j}\,dy_{p,i,j}}{\mathrm{area}_{p,i,j}}. $$

The result is the same first-order dissipative LF transport as
`transport_2d_linrood!` (at first order) and `flux_1d`. The rule is
retained as a debugging / oracle scheme — production transport on the
cubed sphere should use the higher-order PPM-based flux (tracked at
dsc-r1i and blocked on the same schema infrastructure).

## Convergence

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The Layer-B walker harness needs (a) cubed_sphere selector dispatch with
panel-connectivity-aware ghost extension (shared follow-up with
<a href="{{< relref "ppm_edge_cubed_sphere" >}}"><code>ppm_edge_cubed_sphere</code></a>
and
<a href="{{< relref "covariant_laplacian_cubed_sphere" >}}"><code>covariant_laplacian_cubed_sphere</code></a>),
(b) a face-stagger / time-varying binding contract for the per-face
Courant <code>$c</code> at <code>u_edge</code> / <code>v_edge</code>
(shared with the Cartesian
<a href="{{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux.json"><code>lax_friedrichs_flux</code></a>
rule and the limiter Layer-B′ TVD harness), and (c) a Layer-B′
MMS-transport sweep (LF flux + cubed-sphere FV divergence + forward Euler
under Williamson-style transport) rather than a bare reconstruction sweep
— the LF flux alone is non-evolutionary. The proposed asymptotic order
is 1, matching the Cartesian core. Until those features land, the
convergence fixture under
<a href="{{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux_cubed_sphere_xi/fixtures/convergence"><code>fixtures/convergence/</code></a>
declares <code>applicable: false</code> and the rendered convergence
plot is suppressed.
</div>

The Layer-A canonical fixture at
[<code>lax_friedrichs_flux_cubed_sphere_xi/fixtures/canonical/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux_cubed_sphere_xi/fixtures/canonical)
pins the c4 (`Nc = 4`) reference output to within `1e-12` relative
tolerance against the imperative `flux_1d` reference, exercising
all four upwind branches (`c < 0`, `c = 0`, `c > 0`, plus the
boundary case `|c| = 1/2`). Cross-panel ghost extension at the
western (`i = 1`) and eastern (`i = Nc + 1`) panel boundaries is
exercised on the c4 grid. The η sibling's Layer-A fixture is
skip-only and points at the ξ fixture — the algebra is symmetric
under coordinate reflection and the Cartesian-core unit tests in
[`test/test_lax_friedrichs_flux_rule.jl`]({{< param repoURL >}}/blob/main/test/test_lax_friedrichs_flux_rule.jl)
pin the algebraic identity that both axis siblings inherit.

## Reference

- Lax, P. D. & Friedrichs, K. O. (1954). "Systems of conservation laws."
  *Comm. Pure Appl. Math.* 7(1), 159–193 — original LF numerical flux.
- Imperative reference:
  [`flux_1d`]({{< param repoURL >}}/blob/main/src/operators/flux_1d.jl)
  ArrayOp on `CubedSphereGrid` (combines LF face flux + FV divergence
  into a single tendency). The cubed-sphere wrapper additionally
  routes through `_get_courant_xi` / `_get_courant_eta` for boundary-aware
  Courant precomputation and `_match_boundary_fluxes_{xi,eta}!` /
  `_match_rotated_boundary_fluxes!` for conservation across panel seams.
- Sibling rules:
  [`lax_friedrichs_flux`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux.json) (Cartesian
  1D core),
  [`transport_2d`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/transport_2d.json)
  (in-panel 2D specialization),
  [`ppm_edge_cubed_sphere`]({{< relref "ppm_edge_cubed_sphere" >}}) (the
  higher-order panel-boundary edge value used by PPM transport).
