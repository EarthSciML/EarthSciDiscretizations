# The per-tracer lateral-inflow gate

Each advected tracer must be able to carry its **own** lateral boundary concentration.
This is the production requirement for a regional air-quality domain: a chemical
mechanism advecting a dozen species needs O₃ flowing in at ~40 ppb *while* NO flows in
at ~0.4 ppt, through the same walls, on the same wind, in the same document.

It is the two-tracer counterpart of `latlon3d_transport_cwc_regional_inflow`, which gates
the same open-wall PPM stack with a single tracer.

## What it exercises

The verified regional patch — lon `[0, 60]`, lat `[20, 60]`, at **12×9×4** — with its
geometry, winds, metric and initial mass profile taken **verbatim** from
`problems/latlon3d_transport_cwc_regional_inflow.esm`. Halo independence is therefore the
only new thing under test. Two tracers ride the one air-mass flux through the one pair of
inflow rules (`ppm_flux_D_lon_mono_inflow_bc`, `ppm_flux_D_lat_mono_inflow_bc`), each
binding its own halo to the rules' `qbc_{w,e,s,n}` **params** at its own call site:

```
D(Mx*qA, qbcA_w, qbcA_e, wrt: lon)      A: halo 1.0 on all four walls
D(Mx*qB, qbcB_w, qbcB_e, wrt: lon)      B: halo 2.0 on the WEST, 1.0 elsewhere
```

## Why both assertions are load-bearing

They fail in opposite directions, and neither alone is sufficient.

| assertion | measured | what its failure would mean |
|---|---|---|
| `devA` Linf = 0.0, tol 0 (`abs` 1e-12) | **0.0 exactly** | carrying the halo as an operand perturbed the scheme — free-stream/CWC no longer holds bitwise |
| `devB` Linf = 1.8742645277859651, tol 1e-9 | **1.8742645277859651** | B's halo is not read, or not read independently of A's |

Tracer A's halo matches its uniform interior `q == 1`, so every CW84 correction is a
difference of equal values — exactly `0.0` in IEEE — each parabola is the constant 1,
every wall flux collapses to `F = M` bitwise, and `devA = mqA - m` stays exactly 0.

Tracer B's west halo (2.0) does **not** match its interior, and the west wind blows *into*
the domain, so that halo is donated through the wall and `devB` leaves zero.

**`devB` is the non-vacuity half.** `devA` alone is satisfied by a rule that ignores its
halo operands entirely — and, more importantly, by the *old* contract. When `qbc_w` was a
free name read out of the consumer's scope, one model scope held exactly one `qbc_w`; both
tracers necessarily resolved it, and `devB` would sit at `0.0` alongside `devA`. **This
document is unexpressible under the free-name contract**, which is precisely what makes it
the regression gate for the halo being a rule param.

## Relationship to the other regional-inflow cases

- `latlon3d_transport_cwc_regional_inflow` — one tracer, all halos matching: the CWC/free-stream
  invariant itself.
- `problems/_falsifier_regional_inflow.esm` — all four halos perturbed ×1.000001; must
  **fail**, establishing non-vacuity for the CWC case (whose reference is the constant zero).
- **this case** — two tracers, halos deliberately *disagreeing*: gates that a halo is
  bound per call site rather than shared across the model.
