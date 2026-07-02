---
title: "Conformance"
description: "The ast/simulation/convergence categories, the single-pathway doctrine, and how the runners work across the five ESS bindings."
weight: 4
---

The library's value is that one `.esm` rule behaves identically under every
ESS binding (Julia, Python, Rust, TypeScript, Go). The conformance suite under
`tests/conformance/` is what makes that a checked property rather than a hope.

## Single pathway (the doctrine)

This repo contains **zero evaluators, zero rule engines, zero numeric
kernels**. Every number in every golden was produced by an *official ESS
binding runner* through the one canonical pipeline:

```
.esm → parse → import/metaparameter resolution → §9.6.3 rewrite fixpoint → official runner
```

`scripts/runners/*` are thin argument-marshalling wrappers over ESS binding
entry points — they contain no logic of their own. If a binding lacks a
capability (say, a canonical-emit CLI), the fix is an upstream issue in
EarthSciSerialization, **never** a local shim, shadow evaluator, or per-rule
special case here. The archive's retired shadow evaluators are the cautionary
tale: the moment a second pathway exists, "conformance" tests the shim, not
the binding.

Tools that *read* library files — the validator, the docs generator — parse
`.esm` JSON structurally and may call official ESS display/pretty-print APIs,
but never numerically evaluate library math.

## The categories

- **`ast/`** — determinism of lowering. Each case is a small consuming model
  with fixed metaparameter bindings; the golden is the canonical emit of the
  post-lowering AST. All five bindings must reproduce it **byte-identically**
  (esm-spec §9.6.8: two bindings expanding the same file must produce
  byte-identical post-lowering ASTs, or the same rejection). This category
  needs no simulator, so it is the first gate a new binding port clears.

- **`simulation/`** — correctness. Each case routes a problem's inline
  §6.6.5 MMS tests through each binding's official simulation pathway at the
  default resolution. Assertions and tolerances live in the problem file;
  the manifest only routes the runners and pins integrator settings per
  binding (`Tsit5` / `RK45` / `Erk` at tight tolerances, so time-integration
  error stays far below the spatial error being asserted).

- **`convergence/`** — order of accuracy. A manifest per case sweeps
  resolutions via loader-API metaparameter bindings; the committed
  `golden/errors.json` holds the Julia-computed error norms, and
  `scripts/check_convergence_order.py` asserts the observed order from that
  file alone (see [MMS and convergence](/guide/mms-and-convergence/)).

- **`regridding/` and `reprojection/`** — the cross-grid entries' own gates:
  conservation and partition-of-unity as exact invariants, per-pair
  areas/weights and projected points toleranced against goldens.

Manifests declare scope honestly: `reference_binding` (always `julia`),
`bindings_required`, and `scope_excluded` with a reason per excluded binding
(e.g. a rewrite-only port with no simulator). A case whose golden has not been
generated yet carries `"status": "pending-golden"` (or `pending-runners`) —
visible in CI as SKIP, and on the docs pages as a pending callout.

## How the runners work

`./scripts/test-conformance.sh` orchestrates: for each category and case, it
invokes each required binding's runner with the manifest, collects results,
and compares against the goldens. Golden *regeneration* is a separate,
deliberate act — `scripts/regenerate-goldens.sh [category …]`, Julia-only —
and a PR that changes goldens must say *why* (spec change, rule change, or
bug; never "refreshed to green").

## CI

`.github/workflows/conformance.yml` runs the always-on **validate** job with
no bindings installed: JSON-schema validation of every library file, the
policy lint (tags contract L001–L008, with fixtures under
`tests/invalid/lint/`), the lint-fixture expectations, the convergence-order
check, and spellcheck. The per-binding matrix and the cross-binding compare
job attach as the runner work packages land. The docs workflow
(`.github/workflows/docs.yml`) is independent: it regenerates this site's
catalog pages from the library files on every change.
