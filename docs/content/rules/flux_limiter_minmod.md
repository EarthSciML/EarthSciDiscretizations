---
title: "flux_limiter_minmod"
slug: "flux_limiter_minmod"
family: "limiter"
grid_family: "cartesian"
rule_kind: "limiter"
accuracy: "O(dx²) in smooth monotone regions; O(dx) at extrema"
applies_to: "limit(r)"
rule_path: "discretizations/finite_volume/flux_limiter_minmod.json"
description: "Minmod TVD flux limiter (Roe 1986) — the most diffusive of the symmetric second-order TVD limiters."
tags: ["limiter", "tvd", "minmod", "monotonicity"]
---

## Stencil

<figure class="figure">
  <img src="/plots/rules/flux_limiter_minmod-interface.png"
       alt="Slope-ratio interface stencil for minmod — q_{i-1}, q_i, q_{i+1} feeding face i+1/2">
  <figcaption>Three-point slope-ratio stencil at face <code>i+1/2</code>
  (positive-velocity branch). The limiter does not introduce a stencil of
  its own — it consumes the upwind and downwind slopes already required by
  the underlying high-order reconstruction and emits a scalar
  <code>φ(r) ∈ [0, 1]</code> that scales the high-order correction added to
  the low-order (upwind) flux.</figcaption>
</figure>

## Limiter curve

<figure class="figure">
  <img src="/plots/rules/flux_limiter_minmod-stencil.png"
       alt="Minmod limiter φ(r) and Sweby second-order TVD region">
  <figcaption>Minmod limiter <code>φ(r) = max(0, min(r, 1))</code> overlaid on the
  Sweby second-order TVD region. The curve sits on the lower edge of the
  region and passes through <code>(1, 1)</code> — making it the most
  dissipative of the symmetric second-order TVD limiters.</figcaption>
</figure>

## Discrete operator

The minmod limiter enters the standard MUSCL/TVD flux assembly at each cell
interface. With the positive-velocity branch (advection speed `a > 0`):

**1. Slope ratio at face `i+1/2`:**

$$r_{i+1/2} \;=\; \frac{q_i - q_{i-1}}{q_{i+1} - q_i}$$

**2. Apply minmod:**

$$\varphi(r) \;=\; \max\!\bigl(0,\;\min(r, 1)\bigr)$$

**3. Limited high-resolution interface flux:**

$$F_{i+1/2} \;=\; F^{\text{low}}_{i+1/2}
                 \;+\; \varphi(r_{i+1/2})\,
                       \bigl(F^{\text{high}}_{i+1/2} - F^{\text{low}}_{i+1/2}\bigr)$$

where `F_low` is the first-order upwind flux (`a · q_i` here) and `F_high`
is any second-order reconstruction (e.g. Lax–Wendroff, Fromm). The negative-
velocity branch is the mirror image with offsets shifted by one cell.

**4. Equivalent slope-limited form** used by MUSCL reconstruction:

$$\sigma_i \;=\; \varphi(r_i)\,\frac{q_{i+1} - q_i}{\Delta x},
\qquad
q^{L}_{i+1/2} \;=\; q_i + \tfrac{1}{2}\,\Delta x\,\sigma_i$$

so the limiter scales the cell-centered slope used to extrapolate to the
left face state. Both forms produce the same TVD scheme.

### Coefficient layout

| symbol | role | grid location |
|---|---|---|
| `q_{i-1}, q_i, q_{i+1}` | inputs to the slope ratio | cell centers |
| `r_{i+1/2}` | upwind / downwind slope ratio | face |
| `φ(r_{i+1/2})` | scalar limiter ∈ [0, 1] | face |
| `F_{i+1/2}` | limited interface flux (rule output) | face |

### Properties

| Property | Value |
|---|---|
| TVD | yes (Sweby region: `0 ≤ φ(r) ≤ min(2, 2r)`, `φ(1) = 1`) |
| Monotonicity-preserving | yes |
| Symmetric | yes (`φ(r)/r = φ(1/r)`) |
| Smooth-extremum behavior | drops to 1st-order by design |
| Order in smooth monotone regions | 2 |
| Order at extrema / `r ≤ 0` | 1 (the limiter clips the slope to zero) |

Reference: Roe (1986), *Ann. Rev. Fluid Mech.* **18**, eq. (35); Sweby (1984),
*SIAM J. Numer. Anal.* **21**(5), fig. 3. Caller wiring lives in
[`discretizations/finite_volume/README.md`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/README.md).

## Convergence

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The minmod
[<code>fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_minmod/fixtures/convergence)
fixture is currently <code>applicable: false</code> — pointwise
order-of-convergence is the wrong acceptance kind for a TVD limiter, which
intentionally degrades to first order at smooth extrema. The right harness
kind is Layer-B' (monotonicity / TVD norms over initial conditions with
sharp gradients); once that fixture kind lands the page will embed an
empirical plot here. Numeric coverage of the φ(r) checkpoints and TVD on a
smooth + square-wave IC already lives at
[<code>tests/fixtures/flux_limiter_minmod/</code>]({{< param repoURL >}}/blob/main/tests/fixtures/flux_limiter_minmod).
</div>
