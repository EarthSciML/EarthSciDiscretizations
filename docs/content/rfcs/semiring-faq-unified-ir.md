---
title: "A semiring-parameterized FAQ IR for ESS arrayops"
description: "Concrete proposal to generalize the ESS arrayop node into a Functional-Aggregate-Query IR over semirings — unifying tensor contraction (ESM/ESD discretization), relational select-multiply-aggregate (ESI), and the data-dependent index-set construction (mesh topology) that currently must live in imperative grid code."
---

> **Status:** Draft proposal (concrete IR). **Bead:** unassigned.
> **Target repo:** EarthSciSerialization (`packages/EarthSciSerialization.jl`, the
> `arrayop` IR and `tree_walk.jl` evaluator). **Staged in:** EarthSciDiscretizations
> `docs/content/rfcs/` because the authoring session could not push to ESS; relocate
> to `EarthSciSerialization/docs/content/rfcs/` before review.

---

## 1. Summary

ESS already evaluates one specialization of a much broader declarative
operation. This RFC proposes generalizing the `arrayop` node from a
fixed **sum-product contraction over dense, name-matched index sets** into a
**Functional Aggregate Query (FAQ) node parameterized by a semiring**, with three
additive capabilities: data-dependent index sets, value-equality joins, and
content-addressed (Skolem) key construction.

The payoff: a single IR whose specializations are (a) today's tensor/stencil
discretization (ESM/ESD), (b) relational select-multiply-aggregate
(EarthSciInventory's `aggregate(derive(join…))`), and (c) the mesh-topology
construction (edge enumeration, connectivity inversion) that **cannot** be
expressed as an einsum today and is therefore stranded in imperative grid code.
Every existing `arrayop` remains valid as the `sum_product` / dense / no-join
special case — the change is a conservative superset.

## 2. Motivation

Two findings from the ESD unstructured-grid work motivate this.

**(a) ESS's evaluator is already most of the way here.** The MPAS/DUO
nearest-neighbour diffusion rules rely on an ESD-side rewrite
(`_rewrite_unstructured_arrayop!`, `_unstructured_const_arrays`) that flattens a
variable-valence FVM coefficient into a precomputed `coeff` array. Inspection of
`tree_walk.jl` shows this rewrite is **redundant with ESS**: `build_evaluator`
already supports (i) **per-cell dynamic reduction bounds** — expression-valued
contracted ranges like `index(n_edges_on_cell, i) - 1`, expanded per output cell
via `_expand_int_range_dyn` — and (ii) **nested gather** — `index(A, index(B,i,k))`
into both state arrays and `const_arrays`, with arithmetic, resolved by
`_resolve_indices`/`_eval_const_int`. So the FVM coefficient could be evaluated
symbolically from primitive mesh arrays with no flattening at all. The IR, not the
evaluator, is the limiting factor: there is no way to *declare* the semiring, the
data-dependent index set, or the joins that a general unstructured operator needs.

**(b) The relational half already exists next door.** EarthSciInventory (ESI) is a
closed relational-algebra AST (`select, filter, map_dim, derive, join, aggregate,
union`) over categorical index sets, reusing ESM's scalar Expression AST for
per-row math. Its core motif — `Emission = Σ Activity × BaseRate × ∏ adjustments`
= `aggregate(derive(join…))` — is **the same formal object** as the FVM stencil
`Σ_k coeff·(u[nbr]−u[i])`: aggregate, in a semiring, a product of factors each
defined on a subset of the index variables. ESI is that operation over a
*categorical* index space; ESD's discretization is it over a *geometric* one. They
were built as separate formats meeting "at the emissions socket." From the IR's
view the socket is just a change of index space.

This RFC names the common parent and proposes ESS adopt it as the `arrayop` IR, so
ESM/ESD/ESI specialize one evaluator instead of three.

## 3. Background — what `arrayop` is today

The current node (as produced by `discretize.jl` and consumed by `tree_walk.jl`):

```json
{
  "op": "arrayop",
  "reduce": "+",
  "output_idx": ["i"],
  "ranges": { "i": [1, 64], "k": [1, 5] },
  "expr":  { "op": "*", "args": [ /* coeff */, /* operand */ ] },
  "args":  [ "..." ]
}
```

Semantics: for each output tuple in `output_idx`, reduce `expr` over the
*contracted* indices (range keys not in `output_idx`) using `reduce`. Contraction
is **positional**: indices combine by sharing a name. Supported reducers today are
`{+, *, max, min}`. Range bounds may be constants or scalar expressions
(`[1, index(valence, i)]`). Factors are state arrays and injected `const_arrays`
(coords, Fornberg weights, mesh metrics). This is a **dense, name-matched,
sum-product (and a few sibling) contraction** — i.e. einsum plus a handful of
reducers and an escape hatch for dynamic bounds.

## 4. The formalism

A **Functional Aggregate Query** computes, over an index space, an aggregate (in a
commutative semiring `(⊕, ⊗)`) of a product of *factors*, each a function of a
subset of the index variables:

```
out[free] = ⊕_{bound}  ⊗_f  factor_f(vars_f)
```

Specializations:

| Semiring `(⊕, ⊗)` | Factors | Index space | = |
|---|---|---|---|
| `(+, ×)` over ℝ | dense arrays | rectilinear/sparse grid | **einsum / ESD discretization** |
| `(+, ×)` / `(min,+)` / … over ℝ | keyed tables | categorical | **ESI inventory pipeline** |
| `(∨, ∧)` over 𝔹 | relations | tuples | **relational join / existence / dedup** |

The Boolean-semiring row is exactly what einsum lacks and what mesh-topology
construction needs (a join/dedup *is* a Boolean FAQ). The one operation outside
even FAQ is **value invention** — minting a fresh key for each distinct tuple
(numbering the discovered edge set). The declarative answer is a **Skolem
function**: name a thing by its content (`edge(min(u,v), max(u,v))`) rather than by
an allocated counter; ESI already does the bounded form of this with its
`composite` index sets and `pack` expression. Optional dense renumbering is then a
separate `rank` pass, not part of the logic.

## 5. Proposed IR

Generalize `arrayop` with **optional, additive** fields. Absence of every new
field reproduces today's semantics exactly.

```json
{
  "op": "arrayop",
  "semiring": "sum_product",          // NEW: named (⊕, ⊗). Default = today.
  "output_idx": ["i"],
  "ranges": {
    "i": { "from": "cells" },          // NEW: index set by name (or [lo,hi] as today)
    "k": { "from": "edges_of_cell", "of": ["i"] }   // data-dependent / ragged
  },
  "join":  [ { "on": [["e", "edge_id"]] } ],   // NEW: value-equality across factors
  "reduce": "+",                        // retained; = semiring ⊕ when `semiring` absent
  "expr":  { "op": "*", "args": [ /* … */ ] },
  "args":  [ "..." ]
}
```

### 5.1 Semiring
A named `(⊕, ⊗)` pair. Initial registry: `sum_product` (default, = today),
`max_product`, `min_sum` (tropical), `bool_and_or` (relational), plus the
non-semiring statistical reducers ESI needs (`count`, `mean`, `weighted_mean`)
admitted as recognized aggregate forms. `reduce` stays as a shorthand for `⊕` when
`semiring` is omitted, preserving existing files.

### 5.2 Index sets (`from`)
A range value may be a dense interval `[lo, hi]` (today), **or** a reference to a
declared index set. An index set is one of:
- a **dense interval** (grid axis) — ESM `domain.spatial`;
- a **categorical enumeration** — ESI `index_sets` (`county`, `fuelType`);
- a **data-derived set** — materialized from another FAQ (§5.5), e.g. the edge set
  discovered from the face→vertex relation.

`of: ["i"]` expresses a **ragged / dependent** inner set (the edges of cell `i`),
which the evaluator already handles via per-cell dynamic bounds. This single
mechanism unifies ESM grid dims and ESI categorical dims.

### 5.3 Value-equality joins (`on`)
Today factors combine only by sharing an index *name* (positional). `join.on`
adds combination by **value equality of key columns**, subsuming ESI `join` and
making connectivity gathers first-class instead of nested-`index` tricks. A
positional einsum is the degenerate case (join on the shared index itself).

### 5.4 Keyed factors
Unify "const array", "state array", and "ESI table" as one concept: a **keyed
factor** mapping an index tuple to a value. `const_arrays` (mesh metrics, coords)
and ESI `tables` become the same kind of `arg`. No evaluator change for the dense
case; tables are factors keyed by categorical tuples.

### 5.5 Skolem keys and `distinct`
Two primitives close the value-invention gap:
- `{"op": "skolem", "args": ["edge", v_lo, v_hi]}` — a deterministic,
  content-addressed key. Generalizes ESI `pack`.
- `distinct: true` on an index-set-producing `arrayop` under the `bool_and_or`
  semiring — set semantics (dedup) materializing a **data-derived index set**.

Together: enumerate the unique edges (`distinct` Boolean FAQ over faces), name each
by a Skolem key, and expose the result as an index set that a geometric FAQ
(§5.2) consumes. An optional `{"op": "rank"}` assigns dense integers for the
array backend.

## 6. Evaluator changes

Most of this already exists; the deltas are bounded.

| Capability | Status in `tree_walk.jl` | Delta |
|---|---|---|
| Dynamic / ragged reduction bounds | **Present** (`_expand_int_range_dyn`) | none |
| Nested gather into arrays + const_arrays | **Present** (`_resolve_indices`, `_eval_const_int`) | none |
| Reducers `+, *, max, min` | **Present** (`_NK_CONTRACTION`) | parameterize by `semiring` |
| Transcendental ops in bodies (`acos`, `sqrt`, …) | **Present** (`_eval_node`) | none |
| Value-equality joins | absent | resolve `join.on` at build time → gather/merge |
| Named / data-derived index sets | partial (dense + dynamic bound) | index-set registry + materialization |
| Skolem keys / `distinct` / `rank` | absent | new resolve passes (build-time) |

Crucially, the existing model — **build-time unroll → compiled `_Node` tree, with
constants inlined as literals** — is preserved. Joins, Skolem keys, and
data-derived sets are all resolved at build time, producing the same compiled
artifact. Performance characteristics for existing rules are unchanged.

## 7. Worked examples

### 7.1 Today's FVM diffusion (unchanged)
`semiring: sum_product` (default), `reduce: "+"`, ragged `k` bound. Identical to
the current `nn_diffusion_*` arrayop — and, per §2(a), the ESD coefficient flatten
becomes unnecessary because the gathered weight evaluates symbolically.

### 7.2 ESI-style `aggregate(derive(join…))`
```json
{ "op": "arrayop", "semiring": "sum_product",
  "output_idx": ["county", "pollutant"],
  "ranges": { "county": {"from": "county"}, "pollutant": {"from": "pollutant"},
              "src": {"from": "sourceType"}, "fuel": {"from": "fuelType"} },
  "join": [ {"on": [["src","sourceType"],["fuel","fuelType"]]} ],
  "reduce": "+",
  "expr": { "op": "*", "args": [ "activity", "base_rate", "temp_adj" ] },
  "args": [ "activity", "base_rate", "temp_adj" ] }
```
The MOVES running-exhaust contraction as one FAQ over categorical index sets —
ESI expressed in the ESS IR, no new evaluator concepts.

### 7.3 Mesh-edge enumeration (the operation einsum can't do)
```json
{ "op": "arrayop", "semiring": "bool_and_or", "distinct": true,
  "output_idx": ["edge"],
  "ranges": { "f": {"from": "faces"}, "a": {"from": "face_vertices", "of": ["f"]},
              "b": {"from": "face_vertices", "of": ["f"]} },
  "filter": { "op": "<", "args": ["a", "b"] },
  "key":    { "op": "skolem", "args": ["edge", "a", "b"] },
  "expr":   { "op": "true" }, "args": [ "faces" ] }
```
Produces the deduplicated edge index set, content-addressed by vertex pair. A
follow-on `rank` densifies it; `edges_of_cell` (§5.2) is the inversion join. The
geometric FAQ in §7.1 then consumes `edges` as a primitive index set — and the
DUO `area_eff = ¼Σ dc·dv` becomes an ordinary `sum_product` FAQ over it rather
than imperative Julia.

## 8. Schema deltas

Additive only (Draft 2020-12). On the `arrayop` object:
- `semiring`: `string` (enum of registered names). Optional; default `sum_product`.
- `ranges[*]`: allow `{ "from": string, "of"?: string[] }` **in addition to**
  the existing `[lo, hi]` tuple.
- `join`: optional array of `{ "on": [[string, string], …] }`.
- `distinct`: optional `boolean`.
- `key`: optional Expression (Skolem term) for index-set-producing arrayops.
- `filter`: optional Expression predicate (already meaningful for ESI parity).

New Expression ops: `skolem` (variadic), `rank` (unary over an index set), `true`.
Index sets gain a registry entry mirroring ESM `domain` dims and ESI `index_sets`.

## 9. Backward compatibility & migration

- **Strict superset.** Every current `arrayop` is the `sum_product` / dense /
  no-join / no-key case. Files without the new fields are unaffected; the schema
  changes are additive; the evaluator's existing paths are untouched.
- **Staged rollout:** (1) land `semiring` parameter + index-set registry (pure
  refactor of existing reducers); (2) add `join.on` resolution; (3) add `distinct`
  + `skolem` + `rank`. ESD can drop `_rewrite_unstructured_arrayop!` once (1)–(2)
  land and rules reference mesh primitives directly (see §2a).
- **Cross-format:** this is the concrete shape of the "future `earthsci-core`
  shared AST" ESI's spec anticipates — ESM/ESD/ESI would import one IR + one
  evaluator.

## 10. Non-goals, risks, open questions

- **Non-goal:** a query optimizer. FAQ admits worst-case-optimal join planning;
  this RFC proposes the *IR*, not a planner. Build-time unroll is retained.
- **Risk — build-time blowup:** data-derived sets over large meshes unrolled
  eagerly could be slow to *compile*. Mitigation: materialize index sets once and
  cache; consider a streaming (Finch/indexed-streams-style) backend later.
- **Risk — surface bikeshedding:** the JSON keys above are illustrative; the
  semantic generalization is the proposal, not the spelling.
- **Open:** registry and identities of admissible semirings (esp. `mean` /
  `weighted_mean`, which are not semirings — admit as recognized aggregate forms?).
- **Open:** does `distinct`/value-invention belong in the evaluator or in a
  separate mesh/topology preprocessing pass that merely *emits* an ESS index set?
  (The pragmatic factoring most of the field uses.)

## 11. References

- Abo Khamis, Ngo, Rudra. *FAQ: Questions Asked Frequently* (PODS 2016) — the
  semiring-FAQ framework that subsumes einsum and relational aggregation.
- Kovach et al. *Indexed Streams* (POPL 2023) — a fused IR for sparse contraction
  **and** relational joins.
- Finch.jl / TACO — structured/ragged-array compilers generalizing einsum to
  data-dependent index sets.
- Nested Relational Calculus / Datalog± (value invention, the chase) — the
  formal home of Skolem key construction.
- In-repo: `EarthSciSerialization/src/tree_walk.jl` (`build_evaluator`,
  `_resolve_indices`, `_expand_int_range_dyn`); ESS RFC
  `per-cell-metric-binding-eval` (per-cell metric binding); ESD
  `discretizations/finite_difference/nn_diffusion_{mpas,duo}.json` and
  `src/ode_problem.jl` (`_unstructured_*`, the machinery this would retire);
  EarthSciInventory `esi-spec.md` (the relational AST and `pack` Skolem precedent).
