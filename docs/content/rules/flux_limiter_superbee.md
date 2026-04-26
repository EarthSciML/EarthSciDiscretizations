---
title: "flux_limiter_superbee"
slug: "flux_limiter_superbee"
families: "limiter"
grid_families: "cartesian"
rule_kinds: "limiter"
accuracy: "O(dx²) in smooth monotone regions; compressive near discontinuities"
applies_to: "limit(r)"
rule_path: "discretizations/finite_volume/flux_limiter_superbee.json"
description: "Superbee TVD flux limiter (Roe 1986) — sits on the upper edge of the second-order TVD region."
tags: ["limiter", "tvd", "superbee", "compressive"]
---

## Limiter curve

<figure class="figure">
  <img src="/plots/rules/flux_limiter_superbee-stencil.png"
       alt="Superbee limiter φ(r) and Sweby second-order TVD region">
  <figcaption>Superbee φ(r) overlaid on the Sweby (1984) second-order TVD
  region. The curve hugs the <em>upper</em> edge — at r = 1 it is twice the
  identity (slope 2 from the origin), then saturates at 2 — which is what
  makes the scheme compressive: it amplifies sharp gradients but degrades
  smooth extrema more aggressively than minmod.</figcaption>
</figure>

The rule is a scalar AST in the slope-ratio variable `$r`; it carries no
spatial stencil of its own. The caller computes `r` at each interface and
multiplies the high-order slope correction by φ(`r`) — see the [discrete
form](#discrete-form) below and the worked example in
[`discretizations/finite_volume/README.md`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/README.md#composing-a-limiter-with-a-reconstruction).

## Discrete form

Slope ratio at cell interface `i + ½` (1D, uniform Cartesian, `u > 0`):

$$r_i \;=\; \frac{q_i - q_{i-1}}{q_{i+1} - q_i + \varepsilon}$$

(small `ε` guards the locally-flat case). The **superbee** limiter is

$$\varphi(r) \;=\; \max\!\bigl(0,\;\min(2r,\,1),\;\min(r,\,2)\bigr).$$

A MUSCL-style high-order flux scales the slope correction by φ:

$$F_{i+\frac{1}{2}} \;=\; u\,\Bigl(q_i \;+\; \tfrac{1}{2}\,\varphi(r_i)\,(q_{i+1} - q_i)\Bigr).$$

Under forward Euler with `CFL ≤ 1 / (1 + 0.5·φ_max)` (φ_max = 2 for
superbee, so CFL ≤ 0.5), this scheme is strictly TVD. The composition is
reconstruction-agnostic: the same φ(`r`) factor applies when pairing with
PPM or WENO-5 reconstructions.

| Property | Value |
|---|---|
| TVD | yes |
| Monotonicity-preserving | yes — φ(r) = 0 for r ≤ 0 |
| Sweby upper bound | φ(r) ≤ 2 |
| Consistency | φ(1) = 1 |
| Symmetric | yes — φ(r)/r = φ(1/r) |
| Smooth-extremum behavior | compressive (steepens) |
| φ\_max | 2 |

References: Roe (1986), *Ann. Rev. Fluid Mech.* 18:337–365, eq. (36);
Sweby (1984), *SIAM J. Numer. Anal.* 21(5):995–1011, fig. 4.

## Convergence

<div class="callout callout-pending">
<strong>Convergence plot pending fixture activation.</strong>
The convergence fixture at
[<code>discretizations/finite_volume/flux_limiter_superbee/fixtures/convergence/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_superbee/fixtures/convergence)
currently declares <code>applicable: false</code>: a slope-ratio limiter is a
scalar AST, not a stencil, and its acceptance criterion is
<em>monotonicity preservation</em> (Sweby region, φ(1) = 1, strict TVD on
slope-ratio inputs), not asymptotic convergence-order on a manufactured
solution. The Layer-B′ monotonicity / TVD harness lands the matching
fixture kind; until then the
[<code>fixtures/monotonicity/</code>]({{< param repoURL >}}/blob/main/discretizations/finite_volume/flux_limiter_superbee/fixtures/monotonicity)
sweep — Sweby-region properties verified directly off the rule's AST plus
a CFL-bounded TVD check on a smooth + square-wave initial condition — is
the authoritative numeric coverage for this rule.
</div>

## See also

- [`flux_limiter_minmod`]({{< ref "/rules/flux_limiter_minmod" >}}) — the
  diffusive lower-edge counterpart.
- [`ppm_reconstruction`]({{< ref "/rules/ppm_reconstruction" >}}) and
  [`weno5_advection`]({{< ref "/rules/weno5_advection" >}}) — high-order
  reconstructions commonly composed with a TVD limiter.
- [`discretizations/finite_volume/README.md`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/README.md) —
  worked composition example.
