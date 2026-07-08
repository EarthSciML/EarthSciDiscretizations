---
title: "vertical_remap"
slug: "vertical_remap"
families: "finite_volume"
grid_families: "vertical"
rule_kinds: "scheme"
accuracy: "O(dp^3) interior; O(dp) within two cells of the column top/bottom and at active CW84 limiter constraints"
applies_to: "vertical_remap(q, dp_old, dp_new), dim=k"
rule_path: "discretizations/finite_volume/vertical_remap.json"
description: "Conservative PPM vertical remap (Lin 2004) of a cell-averaged scalar from one set of vertical layer thicknesses to another: 4th-order CW84 edge interpolation + Colella–Woodward (1984) §4 monotonicity limiter + exact integration of the sub-grid parabola over the cumulative-pressure overlap. Documentation-only — deferred to a future ESS phase-hook RFC."
tags: ["finite-volume", "ppm", "vertical-remap", "lin-2004", "colella-woodward", "conservative", "monotone"]
---

<div class="callout">
<strong>Documentation-only rule.</strong>
Conservative remap is structurally a <em>phase-hook</em> operation
(Lagrangian → Eulerian re-gridding run between timesteps), not a §7
stencil rule: it has a variable-arity, data-dependent neighbor list (the
cumulative-pressure intersection of <code>dp_old</code> with
<code>dp_new</code>) and binds time-varying per-column metrics that the §7
stencil schema does not encode. This file captures the math, references,
and closed-form AST fragments so the scheme is reproducible, but it is
<strong>not executed by the rule engine</strong> and is excluded from the
walker's Layer-B / Layer-D membership — disposition
<code>deferred-to-phase-hook-rfc</code>. No imperative implementation
exists in any binding (Julia / Python / Rust / TypeScript) and none is to
be re-added. Models needing a Lagrangian vertical configuration must use
an Eulerian-vertical scheme today (e.g.
[`centered_2nd_uniform_vertical`]({{< ref "/rules/centered_2nd_uniform_vertical" >}}))
or wait for the phase-hook RFC.
</div>

## Continuous form

`vertical_remap` conservatively re-maps a cell-averaged scalar $q$ from
an old set of vertical layer thicknesses $dp_{\text{old}}$ to a new set
$dp_{\text{new}}$ (Lin 2004, the "vertically Lagrangian" FV framework).
Working in a pressure-like vertical coordinate, the cumulative interface
pressures are

$$p_{\text{old}}[k+1] = p_{\text{old}}[k] + dp_{\text{old}}[k],\qquad
  p_{\text{new}}[k+1] = p_{\text{new}}[k] + dp_{\text{new}}[k],$$

with reference $p_{\text{old}}[0] = p_{\text{new}}[0] = 0$ (only
differences matter under conservation). Each old layer carries a sub-grid
parabola; the new layer averages are obtained by **exact integration of
those parabolas over the overlap** of each new interval with the old
intervals. Inputs are column-shaped `[Nk]`:

| Input | Meaning |
|---|---|
| `q` | cell-averaged scalar on the old vertical layers |
| `dp_old` | old layer thicknesses (units cancel in `q_new`) |
| `dp_new` | new layer thicknesses; for column-mass invariance $\sum_k dp_{\text{new}}[k] = \sum_k dp_{\text{old}}[k]$ |

The single output is the remapped column `q_new`.

## Reconstruction (PPM edge values)

Layer-interface values use the **4th-order CW84 interior interpolation**
(eq. 1.6) of cell-averaged $q$ — the same coefficients as
[`ppm_reconstruction`]({{< ref "/rules/ppm_reconstruction" >}}):

$$q_{k+1/2}\;=\;\tfrac{7}{12}\big(q_{k}+q_{k+1}\big)\;-\;\tfrac{1}{12}\big(q_{k-1}+q_{k+2}\big).$$

In stencil terms the interior 4-point support and its coefficients are

| selector | axis | stagger | offset | coeff |
|---|---|---|---:|---|
| `vertical` | `$k` | `cell_center` | −1 | −1/12 |
| `vertical` | `$k` | `cell_center` |  0 |  7/12 |
| `vertical` | `$k` | `cell_center` | +1 |  7/12 |
| `vertical` | `$k` | `cell_center` | +2 | −1/12 |

The accuracy degrades near the column ends as the 4-cell support runs
out:

| Layer | Edge-value rule | Order |
|---|---|---|
| $k = 1$ or $k = N_k$ | $q_{\text{left}}[k] = q_{\text{right}}[k] = q[k]$ | O(1) constant on the boundary layer |
| $k = 2$ or $k = N_k-1$ | $q_{\text{left}}[k] = \tfrac12(q_{k-1}+q_k)$, $q_{\text{right}}[k] = \tfrac12(q_k+q_{k+1})$ | O(dp²) centered average, then CW84 limiter |
| interior | 4-point CW84 stencil above, then CW84 limiter | O(dp⁴) edge value, O(dp³) reconstruction |

## Monotonicity limiter

Each layer's $(q_{\text{left}}, q_{\text{right}}, q)$ triple is passed
through the **Colella & Woodward (1984) §4 monotonicity limiter**,
encoded as a closed-form `ifelse` AST so the rule is evaluator-agnostic.
This is the **same limiter block** used by
[`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}}). With
$dq = q_R - q_L$ and $q_6 = 6\big(q - \tfrac12(q_L + q_R)\big)$:

- **`is_extremum`** — when $(q_R - q)(q - q_L) \le 0$, $q$ is a local
  extremum within the column; flatten the parabola to the constant $q$.
- **`overshoot_left`** — when $dq\cdot q_6 > dq^2$ the right edge sits
  past the parabola's interior extremum; pull the **left** edge inward to
  $q_L = 3q - 2q_R$.
- **`overshoot_right`** — when $-dq^2 > dq\cdot q_6$, pull the **right**
  edge inward to $q_R = 3q - 2q_L$.

## Conservative remap sweep

Within old layer $k$ the limited $(q_L, q_R, q)$ define the CW84 sub-grid
parabola (eqs. 1.5, 1.7, 1.10) in a local coordinate $\xi \in [0,1]$
running from the lower interface ($\xi = 0$) to the upper ($\xi = 1$):

$$a(\xi)\;=\;a_L \;+\; \xi\big(\,da + a_6(1-\xi)\,\big),\qquad
  da = a_R - a_L,\quad a_6 = 6\big(q - \tfrac12(a_L + a_R)\big).$$

Its closed-form indefinite integral (Lin 2004 eq. 16) is

$$F(\xi)\;=\;a_L\,\xi \;+\; \tfrac12\,\xi^2\,(da + a_6)\;-\;\tfrac{a_6}{3}\,\xi^3.$$

The remap sweep is **not a fixed stencil** over neighbors: for each new
layer $k_n$, integrate the parabolas of every old layer $k_o$ whose
interface span $[p_{\text{old}}[k_o-1], p_{\text{old}}[k_o]]$ intersects
$[p_{\text{new}}[k_n-1], p_{\text{new}}[k_n]]$, with local coordinate
$\xi(p, k_o) = (p - p_{\text{old}}[k_o-1])/dp_{\text{old}}[k_o]$:

$$\text{integral}[k_n] = \sum_{k_o:\,\text{overlap}}
  \Big(F\big(\xi(\min(p_{\text{new}}[k_n],\,p_{\text{old}}[k_o]),\,k_o)\big)
     - F\big(\xi(\max(p_{\text{new}}[k_n-1],\,p_{\text{old}}[k_o-1]),\,k_o)\big)\Big)\,dp_{\text{old}}[k_o],$$

$$q_{\text{new}}[k_n] = \text{integral}[k_n] / dp_{\text{new}}[k_n].$$

A degenerate $dp_{\text{new}}[k_n] \to 0$ falls back to
$q_{\text{new}}[k_n] = q[k_o]$ for the bracketing old layer.

### Conservation invariant

Column-integrated mass is preserved exactly in real arithmetic
(floating-point error bounded by $O(N_k\,\varepsilon\sum_k |q_k|\,dp_{\text{old}}[k])$):

$$\sum_{k_n} q_{\text{new}}[k_n]\,dp_{\text{new}}[k_n]
  \;=\; \sum_{k_o} q[k_o]\,dp_{\text{old}}[k_o]
  \qquad\text{when}\quad \textstyle\sum_{k_n} dp_{\text{new}}[k_n] = \sum_{k_o} dp_{\text{old}}[k_o].$$

## Convergence

As a documentation-only phase-hook scheme, `vertical_remap` is excluded
from the walker's Layer-B membership and ships no MMS convergence sweep.
Its convergence fixture under
[`discretizations/finite_volume/vertical_remap/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/vertical_remap/fixtures/convergence)
documents the intended accuracy (O(dp³) interior) but is suppressed
pending the phase-hook RFC and its associated remap runner.

## Reference

- Lin, S.-J. (2004). "A 'Vertically Lagrangian' Finite-Volume Dynamical
  Core for Global Models." *Monthly Weather Review*, 132(10), 2293–2307,
  eqs. (16)–(18) — conservative remap framework.
- Colella, P. and P. R. Woodward (1984). *Journal of Computational
  Physics*, 54(1), 174–201, eqs. (1.6)–(1.10) — PPM reconstruction and
  the (1.7)–(1.10) monotonicity limiter.
- Shares its edge interpolation with
  [`ppm_reconstruction`]({{< ref "/rules/ppm_reconstruction" >}}) and its
  monotonicity limiter with
  [`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}}).
