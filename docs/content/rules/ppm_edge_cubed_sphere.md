---
title: "ppm_edge_cubed_sphere"
slug: "ppm_edge_cubed_sphere"
families: "finite_volume"
grid_families: "cubed_sphere"
rule_kinds: "scheme"
accuracy: "O(h²) at panel-boundary interfaces"
applies_to: "reconstruct_panel_edge(q)"
rule_path: "discretizations/finite_volume/ppm_edge_cubed_sphere.json"
description: "FV3 two-sided PPM extrapolation at cubed-sphere panel-boundary interfaces — averages two one-sided edge values to bridge the discontinuous panel coordinate system."
tags: ["finite-volume", "cubed-sphere", "ppm", "panel-boundary", "fv3"]
---

## Continuous form

At a cube-edge interface the local panel coordinate system is discontinuous:
`(ξ, η)` on one panel does not extend smoothly into `(ξ', η')` on the
neighbor panel. The standard 4th-order PPM edge formula
([CW84]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_reconstruction.json) eq. 1.6)

$$ q_{i+1/2} \;=\; \tfrac{7}{12}(q_i + q_{i+1}) \;-\; \tfrac{1}{12}(q_{i-1} + q_{i+2}) $$

assumes a single uniform coordinate system over the four contributing
cells; FV3 §6.5 replaces it at panel boundaries with a two-sided formula
that respects the cell widths on each side independently:

$$ a_E \;=\; \tfrac{1}{2}\,\frac{(2\,dx_1 + dx_2)\,q_1 - dx_1\,q_2}{dx_1 + dx_2} \;+\; \tfrac{1}{2}\,\frac{(2\,dx_0 + dx_{-1})\,q_0 - dx_0\,q_{-1}}{dx_0 + dx_{-1}}. $$

`q_{-1}, q_0` are the two cells on the "minus" side of the interface
(neighbor panel through ghost cells when the interface is the panel
seam); `q_1, q_2` are the two cells on the "plus" side (local panel).
On the gnomonic *equidistant* cubed sphere the spacing is uniform
isotropic `h = dξ = dη = π/(2 Nc)` and the formula collapses to a
linear 4-cell stencil with rational coefficients.

## Discrete operator

Substituting `dx_{-1} = dx_0 = dx_1 = dx_2 = h` into the FV3 two-sided
formula gives, along axis `$x ∈ {xi, eta}`,

$$ q_{i+1/2} \;=\; -\tfrac{1}{4}\,q_{i-1} \;+\; \tfrac{3}{4}\,q_{i} \;+\; \tfrac{3}{4}\,q_{i+1} \;-\; \tfrac{1}{4}\,q_{i+2}. $$

The 4-cell stencil at offsets `{-1, 0, 1, 2}` produces the value at the
panel-edge interface lying between offsets 0 and 1:

| offset along `$x` | coefficient | role |
|---|---|---|
| `-1` | `-1/4` | minus-side outer (`q_{-1}`) |
| ` 0` | `+3/4` | minus-side inner (`q_{0}`) |
| `+1` | `+3/4` | plus-side inner (`q_{1}`) |
| `+2` | `-1/4` | plus-side outer (`q_{2}`) |

Compare with the standard interior PPM edge formula
([`ppm_reconstruction.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_reconstruction.json),
offsets `{-1, 0, 1, 2}` with coefficients `[-1/12, 7/12, 7/12, -1/12]`):
same 4 cells, different weights. The two-sided formula is exact for
linear fields only — by design, since the higher-order interior
formula's accuracy claim no longer holds across the panel discontinuity.

Cross-panel offsets resolve to neighbor-panel ghost cells (with the
appropriate 2×2 basis rotation) through the cubed_sphere grid accessor
(`src/grids/panel_connectivity.jl`); the rule selectors do **not**
carry a `panel` field per [SELECTOR_KINDS.md decision #13]({{< param repoURL >}}/blob/main/discretizations/SELECTOR_KINDS.md).

### Monotonicity constraint (FV3 eq. 6.6)

The FV3 monotonicity guard clamps the linear stencil result to the
range of the four contributing cells,

$$ a_E := \mathrm{clamp}\bigl(a_E,\; \min(q_{-1}, q_0, q_1, q_2),\; \max(q_{-1}, q_0, q_1, q_2)\bigr), $$

preventing new extrema at panel-boundary interfaces. This is a
**non-linear post-processing step**, not a stencil coefficient — the
rule carries the linear stencil only and documents the clamp under
its `monotonicity_constraint` block. Encoding post-stencil clamps in
the ESS §7 schema (parallel to the AST `formula` block in
[`flux_limiter_minmod.json`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_minmod.json))
is a tracked follow-up.

## Convergence

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The Layer-B walker harness needs the in-flight cubed_sphere dispatch
with panel-connectivity-aware ghost extension to evaluate this rule on
a manufactured solution defined in <code>(ξ, η)</code>. Until that
extension lands, the convergence fixture under
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_edge_cubed_sphere/fixtures/convergence)
declares <code>applicable: false</code> and the rendered convergence
plot is suppressed. The proposed asymptotic order is 2 — the two-sided
formula is exact for linear fields by construction, and the
truncation error on smooth fields is `O(h²)` from the omitted
quadratic term.
</div>

The bit-equivalence canonical fixture at
[<code>fixtures/canonical/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/ppm_edge_cubed_sphere/fixtures/canonical)
pins the c24 (`Nc = 24`) reference output to within `1e-12` relative
tolerance against the pre-port `ppm_edge_value_twosided` imperative
reference (see "Reference" below).

## Reference

- Harris et al. (2021), *GFDL FV3 Technical Memorandum*, Eq. 6.5–6.6 —
  two-sided extrapolation at cube edges and the monotonicity guard.
- Imperative reference: `ppm_edge_value_twosided` /
  `ppm_edge_value_twosided_limited` in
  `src/operators/ppm_edge.jl` (deleted at port time; the bit-equivalence
  canonical fixture is the lasting contract).
