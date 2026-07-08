---
title: "lax_friedrichs_flux"
slug: "lax_friedrichs_flux"
families: "finite_volume"
grid_families: "cartesian"
rule_kinds: "scheme"
accuracy: "O(h)"
applies_to: "flux(q, c), dim=x"
rule_path: "discretizations/finite_volume/lax_friedrichs_flux.json"
description: "Lax–Friedrichs first-order numerical flux for 1D linear advection on a uniform Cartesian axis. For face-staggered Courant c the LF flux reduces to first-order upwinding F_{i+1/2} = max(c,0)·q_i + min(c,0)·q_{i+1}, with upwind cell selection encoded directly in an abs op (no caller-side branching on sign(c))."
tags: ["finite-volume", "lax-friedrichs", "upwind", "flux", "monotone", "tvd", "first-order"]
---

## Continuous form

`lax_friedrichs_flux` discretizes the conservative 1D advective face flux
on a uniform Cartesian grid for linear advection $f(q) = c\,q$,

$$F(x_{i+1/2},\,t)\;=\;f\big(q(x_{i+1/2},\,t)\big),$$

with a **face-staggered Courant number** $c_{i+1/2} = u_{i+1/2}\,\Delta t/\Delta x$.
The Lax & Friedrichs (1954) numerical flux is

$$F_{i+1/2}\;=\;\tfrac{c}{2}\,(q_i + q_{i+1})\;-\;\tfrac{|c|}{2}\,(q_{i+1} - q_i),$$

which, collecting terms, is algebraically identical to first-order
upwinding:

$$F_{i+1/2}\;=\;\frac{c + |c|}{2}\,q_i \;+\; \frac{c - |c|}{2}\,q_{i+1}\;=\;\max(c,0)\,q_i \;+\; \min(c,0)\,q_{i+1}.$$

Upwind cell selection is encoded **directly in the `abs` op** — no
caller-side branching on $\operatorname{sign}(c)$ is required.

## Stencil

A two-point face stencil straddling the interface $F_{i+1/2}$ (the right
cell of the interface is offset `0`):

| selector kind | axis | offset | coeff |
|---|---|---:|---|
| `cartesian` | `$x` | 0 | `(c + |c|) / 2` |
| `cartesian` | `$x` | 1 | `(c − |c|) / 2` |

Combined as `+`. For $c > 0$ the offset-`0` coefficient is $c$ and the
offset-`1` coefficient is $0$ (pure upwind from the left cell); for
$c < 0$ the weights swap (upwind from the right cell). The closed
`flux_form` AST in the JSON states the same identity directly as

$$F_{i+1/2}\;=\;\frac{c + |c|}{2}\,q_i \;+\; \frac{c - |c|}{2}\,q_{i+1}.$$

## Properties

- **Monotone** and **TVD** — the scheme is dissipative by construction; it
  adds the first-order numerical diffusion characteristic of upwinding and
  cannot produce new extrema.
- The face-Courant `$c` is a **per-face binding** supplied by the caller
  (analogous to the slope-ratio `$r` in
  [`flux_limiter_minmod`]({{< ref "/rules/flux_limiter_minmod" >}})), not
  a stencil-offset selector.
- Retained as a **debugging / oracle scheme**. Production transport should
  use the higher-order PPM or WENO rules
  ([`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}}),
  [`weno5_advection`]({{< ref "/rules/weno5_advection" >}})).

## Composition

The caller computes the face-staggered Courant
$c_{i+1/2} = u_{i+1/2}\,\Delta t/\Delta x$ and applies this rule at every
interface; the cell-centred tendency follows from the FV divergence

$$T_i\;=\;-\,\frac{F_{i+1/2} - F_{i-1/2}}{h}$$

on a uniform grid (or composes with
[`divergence_arakawa_c`]({{< ref "/rules/divergence_arakawa_c" >}}) on a
C-grid). It is the **first-order sibling** of
[`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}}): same `op = flux` shape
and per-face binding contract, with the LF `abs`-based first-order upwind
in place of PPM's higher-order reconstruction + Courant-fraction flux
integral.

## Convergence

The rule's Layer-B convergence fixture lives under
[`discretizations/finite_volume/lax_friedrichs_flux/fixtures/convergence/`]({{< param repoURL >}}/blob/main/discretizations/finite_volume/lax_friedrichs_flux/fixtures/convergence).
Like [`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}}), a flux-form
transport sweep needs the per-face Courant binding (`$c`) and a Layer-B′
MMS-transport runner (numerical flux + FV divergence + time stepping) to
measure empirical order on a smooth profile; until those harness features
land the convergence fixture is suppressed and the algebraic upwind
identity above is exercised inline.

## Reference

- Lax, P. and K. O. Friedrichs (1954). "Weak solutions of nonlinear
  hyperbolic equations and their numerical computation." *Comm. Pure
  Appl. Math.* 7(1), 159–193 — original LF flux.
- First-order reduction: for linear advection with face-staggered Courant
  $c$, $F_{i+1/2} = \max(c,0)\,q_i + \min(c,0)\,q_{i+1}$.
- Higher-order sibling: [`flux_1d_ppm`]({{< ref "/rules/flux_1d_ppm" >}})
  (PPM flux-form transport).
