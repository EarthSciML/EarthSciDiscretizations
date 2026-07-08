---
title: "ArrayDiscretization parity contract"
description: "Authoritative reference for the MOL PR #531 pin, the parity definition, the deliberate frame-split divergence, and the golden-capture convention for all array-discretization-parity work."
bead: "esd-cm7"
---

## Purpose

This document is the authoritative anchor for all work in the
`array-discretization-parity` label group. Subsequent rule-authoring and
fixture beads cite this document for:

- the frozen MOL PR #531 SHA they compare against,
- the precise definition of parity,
- the deliberate divergence that reviewers must not flag as a parity gap, and
- the convention for capturing and regenerating conformance goldens.

---

## 1. Pinned reference

All `array-discretization-parity` work is pinned to:

| Field | Value |
|---|---|
| Repository | [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) |
| Pull request | [#531](https://github.com/SciML/MethodOfLines.jl/pull/531) |
| **Pinned SHA** | **`35cc9143dc553ac7d3619738bd77b250c1ed162f`** |

**PR #531 is open and may change.** All parity work and golden captures
freeze against the SHA above, not the live branch tip. Do not update this
SHA without a deliberate decision bead.

---

## 2. Parity definition

> **Parity** = the ESS `discretize(esm)` pipeline, fed ESD catalog rules,
> reproduces #531's **interior** ArrayOp output — validated as
> byte-identical canonical AST across bindings
> (ESS `discretization.md` §13), captured **once** from the pinned SHA and
> frozen.

Three properties of this definition are load-bearing:

1. **Interior only.** The comparison scope is the interior ArrayOp — the
   cells away from domain edges. Boundary and periodic handling is
   delegated to ESD's boundary/periodic rules (see §4 on deliberate
   divergence).

2. **Byte-identical canonical AST.** The conformance check is not
   floating-point tolerance on evaluated values; it is structural identity
   of the lowered AST as specified in ESS `discretization.md` §13. Every
   binding must agree on the same canonical form.

3. **No live MOL dependency.** The golden is captured once at the pinned
   SHA and committed to `tests/conformance/`. Subsequent CI runs compare
   against the committed golden; they do not re-run MOL or fetch from
   GitHub. The ESD test suite has no runtime dependency on MOL or any
   external network resource.

---

## 3. No system-level discretizer in ESD

ESD does **not** carry a `discretize` function or rule evaluator analogous
to `MethodOfLines.discretize`. The MOL-`discretize` analogue is
**ESS's `discretize(esm)`**, which is already landed.

ESD's contribution to the pipeline is:

- **catalog rules** — JSON rule files under `discretizations/` that express
  interior stencils as closed `arrayop` ASTs in the ESS §4.2 op vocabulary,
- **conformance fixtures** — pinned goldens under `tests/conformance/` that
  verify cross-binding agreement on those rules.

ESD never grows a parallel evaluator, a shadow rule runner, or a binding-local
fast path that bypasses ESS. Any such addition is an AGENTS.md single-pathway
violation and must be retired, not extended.

---

## 4. Deliberate divergence

ESD/ESS emit **one** `arrayop` over the interior. They do **not** replicate
#531's internal frame / full-interior mode split.

That split is a **MOL scalarizer optimization**: #531 distinguishes a "frame"
(boundary-adjacent cells that straddle the edge of the interior stencil) from
the "full interior" (cells whose entire stencil lies in the interior) and emits
separate lowerings for each. This separation is not required by the mathematical
parity contract — it is an artifact of the MOL lowering strategy.

ESD relies instead on:

- A **single interior rule** that expresses the full stencil as an `arrayop`,
- **Boundary rules** (e.g. `periodic_bc`, `dirichlet_bc`) that rewrite the
  AST at boundary-adjacent cells, where some stencil neighbors fall outside
  the interior.

The net result is mathematically equivalent output over the interior, but
the AST structure differs from #531's frame/full-interior split.

**Reviewers must not flag this as a parity gap.** The parity contract (§2) is
interior-only ArrayOp agreement — not AST isomorphism with #531's internal
lowering stages.

---

## 5. Golden-capture convention

### 5.1 One-time capture at the pinned SHA

Capture goldens **once**, from the pinned SHA `35cc9143dc553ac7d3619738bd77b250c1ed162f`.
Do not re-capture against HEAD of #531; the branch is open and may change.

```bash
# In a local MethodOfLines.jl checkout:
git checkout 35cc9143dc553ac7d3619738bd77b250c1ed162f

# Run the ESD canonical pipeline against the pinned MOL state,
# capturing the interior ArrayOp output.
# (exact invocation defined in the per-fixture README under tests/conformance/)
```

### 5.2 Storage location and header convention

Goldens live under `tests/conformance/`, mirroring the existing conformance
layout:

```
tests/conformance/
  mol531/
    <rule-name>/
      README.md          ← describes the fixture, cites the pinned SHA
      fixtures.json      ← fixture inputs
      golden/
        <fixture>.json   ← captured output
      regenerate_golden.jl   ← (or .py / .mjs per binding)
```

Every fixture file and README **must** include the pinned SHA in its header
so the capture's provenance is self-documenting:

```json
{
  "_mol531_sha": "35cc9143dc553ac7d3619738bd77b250c1ed162f",
  "_captured_by": "esd-<bead-id>",
  ...
}
```

### 5.3 Regeneration via the canonical pipeline

Goldens are regenerated via the canonical ESS pipeline — rule application in
ESS → ArrayOp → eval through the ESD passthrough — **not** via a
binding-local reimplementation. This matches the rule from AGENTS.md and
from the existing `tests/conformance/grids/*/regenerate_golden.*` scripts.

A `regenerate_golden.*` script that computes expected values via a
binding-local code path is an AGENTS.md single-pathway violation; the
resulting golden masks divergence between bindings rather than exposing it.

### 5.4 Fixture scope

Each per-rule golden captures the **interior** ArrayOp output at the input
conditions used by #531's test suite at the pinned SHA. The fixture scope
is intentionally narrow — enough to validate parity, not a full regression
suite.

---

## Cross-references

- **AGENTS.md** — single-pathway rule (ESD is a catalog over ESS, not a
  parallel runtime).
- **ESS `discretization.md` §13** — canonical AST byte-identity protocol
  for cross-binding conformance.
- **`docs/rule-catalog-roadmap.md` §A** — MOL PR #531 source section in
  the Phase-0 catalog inventory (pinned SHA also recorded there).
- **Label `array-discretization-parity`** — all rule and fixture beads in
  this effort; each cites this document for the contract.
