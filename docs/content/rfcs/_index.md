---
title: "RFCs"
description: "Authoritative design notes and parity contracts for ESD development efforts."
---

Short RFC-style notes that establish contracts, pin external references, and
record deliberate design decisions that downstream beads must cite rather than
re-derive.

- **[ArrayDiscretization parity contract]({{< ref "/rfcs/array-discretization-parity" >}})** —
  MOL PR #531 frozen SHA, parity definition, deliberate frame-split divergence,
  and golden-capture convention for the `array-discretization-parity` effort.
- **[A semiring-parameterized FAQ IR for ESS arrayops]({{< ref "/rfcs/semiring-faq-unified-ir" >}})** —
  proposal (targeting EarthSciSerialization) to generalize the `arrayop` node into a
  Functional-Aggregate-Query IR over semirings, unifying ESM/ESD tensor contraction,
  ESI relational select-multiply-aggregate, and data-dependent index-set construction.
  Staged here pending relocation to the ESS repo.
