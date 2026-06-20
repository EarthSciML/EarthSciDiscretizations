---
title: "Implementing the build-time relational engine across ESS bindings"
description: "Research note and implementation recommendation for the setup-time relational engine (distinct / join / skolem / rank / group-by) that the semiring-FAQ IR requires, across the Julia, Rust, and Python implementations of earthsci-toolkit. Prioritizes existing libraries and pins a cross-binding determinism spec."
---

> **Status:** Research note / implementation recommendation. Companion to
> [*A semiring-parameterized FAQ IR for ESS arrayops*](semiring-faq-unified-ir.md)
> (§6.1 static/dynamic partition, §9 v1 topology engine), which introduces the
> need for a build-time relational engine in `earthsci-toolkit`.
> **Target repo:** EarthSciSerialization. **Staged in:** EarthSciDiscretizations
> `docs/content/rfcs/` for the same reason as its parent RFC.

---

## 1. What this engine is, and the constraint that shapes it

The unified IR makes mesh topology first-class: `distinct` (edge enumeration),
`join` (connectivity inversion), `skolem` (content-addressed keys), and `rank`
(dense renumbering) become declarable operations (RFC §5.5). The static/dynamic
partition (RFC §6.1) runs these **once at setup**, off the per-timestep hot path,
to materialize index sets that the numeric stencil then consumes. The engine needs
exactly five primitives over integer-keyed tuples (vertex/cell IDs, scale
10⁴–10⁷, one-time):

1. **`distinct`** — deduplicate tuples (unique mesh edges from face→vertex lists).
2. **`join`** — value-equality equi-join (connectivity inversion, *edges of cell i*).
3. **`skolem`** — deterministic content-addressed key from a tuple.
4. **`rank`** — dense integer renumbering of a distinct set (assign IDs to deduped tuples).
5. **group-by + semiring aggregate** (sum/min/max) over those sets.

**The decisive constraint comes from the ESS architecture.** `earthsci-toolkit` is
**not one core with FFI bindings** — it is *parallel native implementations* per
language (`EarthSciSerialization.jl`, `earthsci-toolkit-rs`, `earthsci_toolkit`
(Python), plus TS/Go), each conforming to a shared contract and verified by a
cross-binding conformance suite (`scripts/test-conformance.sh` against
`CONFORMANCE_SPEC.md`). Therefore the hard problem is **bit-for-bit determinism
across three independent implementations**: identical deduped sets, identical dense
IDs, identical skolem keys. That, not raw speed, governs every choice below.

## 2. Recommendation in one line

**Do not adopt a heavy relational library or a shared embedded engine. In each
binding, hand-roll the five primitives on the data structure that language already
depends on, and enforce a single cross-binding *determinism spec* (canonical
sort-based ordering + canonical-tuple skolem keys).** The relational logic is
~100 lines per language; the value is in the spec, which no library provides for
free.

This also matches what the ESS codebase *already does*: the primitives exist
informally in two of the three bindings and only need to be unified and made
deterministic (§4).

## 3. Why not a shared engine (DuckDB / Polars / Arrow)

A shared engine is the obvious "identical semantics for free" idea. It was
investigated and rejected.

| Engine | Julia | Python | Rust | Same core? | Weight | Verdict |
|---|---|---|---|---|---|---|
| **DuckDB** (1.5.2, 2026) | `DuckDB.jl` ✅ | `duckdb` ✅ | `duckdb-rs` ✅ | **Yes** (one C++ engine via C API) | Heavy: native `libduckdb` ~25–60 MB/platform | **Rejected** — only true 3-language same-core option, but contradicts the parallel-native architecture, adds a heavy native dep to all three packages, and *still requires* `ORDER BY`-everything discipline + conformance tests. Buys less than it appears. |
| **Polars** | `Polars.jl` ✗ (~30★, unmaintained since ~2023) | `polars` ✅ | `polars` crate ✅ | Shared Rust core in **2 of 3** | Medium–heavy | **Rejected** — no mature Julia binding; can't anchor 3-way conformance. (`maintain_order` has also had optimizer escapes.) |
| **Arrow / Acero** | `Arrow.jl` ✗ (format/IO only, no Acero) | `pyarrow` ✅ (Acero) | `arrow-rs` ✗ (uses DataFusion) | **No** (three different execution stacks) | Medium | **Rejected** — not uniform across the three. |

DuckDB remains worth keeping as a **throwaway oracle** during conformance-test
development (`SELECT DISTINCT … ORDER BY …`, `dense_rank() OVER (ORDER BY …)`) to
cross-check the hand-rolled output — but not as a shipped dependency.

## 4. Per-binding recommendation (and what already exists)

Each binding should add a small `relational` module using the structure it
*already* depends on. The cross-binding agreement comes from §5, not the library.

### Julia — stdlib `Dict`/`Set`/`sort` (+ `OrderedCollections`, already a dep)
`EarthSciSerialization.jl` (v0.6.0, Julia ≥1.10) deps are lean — `JSON3`,
`OrderedCollections`, `Tullio`, `Unitful`; **no DataFrames/DuckDB**. The evaluator
is two-phase: build-time `_compile` / `build_evaluator` and hot-path `_eval_node`;
the engine slots in as a build-time topology pre-pass beside `_compile`. **`src/graph.jl`
already hand-rolls distinct (`Set{String}`), joins (composite-key strings), and
group-by (`Dict` node-maps) — but with no `rank` and no ordering guarantees.** The
work is to formalize that and pin the order. `sort` is stable by default since
Julia 1.9. *Reject* DataFrames.jl (multi-second TTFX, "undefined" join/group order)
and DuckDB.jl (native binary) for production.

### Rust — `indexmap` (already a dep) + canonical sort
`earthsci-toolkit-rs` (v0.6.0, edition 2024) already depends on
`indexmap = "2"`, `ndarray`, `smallvec`; **no polars/datafusion/arrow**. It already
encodes the exact discipline needed: `src/performance.rs` does sort-then-enumerate
for reproducible dense indices (**that is `rank`**), and `src/canonicalize.rs`
sorts args on a stable key and has `format_canonical_float` (**that is
`skolem`/`distinct`**). Add `src/relational.rs` using `IndexSet`/`IndexMap` (whose
iteration order is insertion order, *independent of the hasher*) and
`sort_unstable` on the full tuple. *Reject* `polars`/`datafusion`/`arrow` (heavy,
out of proportion) and never let a non-portable fast hasher (`ahash`,
`rustc-hash`/FxHash) drive emitted order or keys.

### Python — NumPy (already a dep) `lexsort`/`unique`/`searchsorted`
`earthsci_toolkit` already hard-depends on NumPy/SciPy/xarray; **no
pandas/polars/duckdb/pyarrow**. The evaluator is a NumPy AST interpreter
(`numpy_interpreter.py`); spatial/mesh ops are contractually lowered at setup
(`UnreachableSpatialOperatorError`) — exactly this engine's slot. Relational code is
**greenfield** here. Build the five primitives on `np.unique(axis=0)` (lexicographically
sorted unique rows), `np.lexsort`, and `searchsorted`-based joins; reuse the
existing `canonicalize.py` total order. *Reject* pandas (dtype coercion, shifting
sort defaults) and bare `set`/`hash()` (PYTHONHASHSEED-sensitive).

## 5. The cross-binding determinism spec (the normative deliverable)

This is the actual product of the research: a spec all three implementations must
honor so their outputs are byte-identical. **Governing principle: every emitted set
is a pure function of a defined total order over tuples; no observable output ever
depends on hash-table iteration order or a language-native hash value.**

1. **Total order.** Lexicographic over tuple fields, documented per type: integers
   by value; **strings by Unicode code-point (UTF-8 byte) order**, not locale
   collation; floats *forbidden in keys* (keep keys integer IDs) — if unavoidable,
   reuse `canonicalize`'s float formatting and reject NaN, normalize `-0.0`→`0.0`.
2. **`distinct`** = sort by the total order, drop adjacent duplicates. Output order
   is the sorted order — **never insertion / first-seen order** (not portable: Rust
   `HashSet` is randomly seeded, Julia `Dict` order is unspecified, Python `set`
   order is PYTHONHASHSEED-sensitive).
3. **`rank`** = dense IDs by position in the sorted distinct sequence. **Pin the
   base in `CONFORMANCE_SPEC.md`** (Julia is 1-based, Rust/Python 0-based): assert
   on the canonical numbering and convert at the binding boundary.
4. **`skolem`** = a **canonical tuple**, not a hash. For an undirected edge use
   `(min(u,v), max(u,v))`; generalize to "sort components for symmetric relations,
   preserve order for directed." The dense ID then comes from `rank`. This keeps
   hashing out of the determinism-critical path entirely.
   - *If* a fixed-width fingerprint is genuinely required, specify a **portable,
     seed-pinned, non-cryptographic hash** (XXH3-64, seed 0 — native impls in all
     three: `xxhash-rust`, `xxhash`, `XXhash.jl`; or BLAKE2/SHA via `hashlib`) over
     a **canonical byte serialization** (fixed field order, little-endian ints,
     UTF-8 strings, length-prefixed). **Never** a language-native `hash()` /
     `Base.hash` / `DefaultHasher`.
5. **`join` / group-by aggregate.** Use hashing only to *bucket*; **sort the
   emitted result by the canonical key**. Semiring combines must be associative +
   commutative (sum, min, max, count, boolean-or) so input/parallel order can't
   change results. Beware floating-point parallel reduction (last-ULP drift): keep
   exact/integer semirings, or do the final reduce sequentially in canonical order.

### The randomization footguns this neutralizes
- **Rust:** `HashMap`/`HashSet` default to SipHash-1-3 with a per-instance random
  seed → "arbitrary" iteration order.
- **Python:** `hash()` of str/bytes is SipHash with per-process `PYTHONHASHSEED`.
- **Julia:** `Dict`/`Set` iteration order is an unspecified implementation detail;
  `Base.hash` is process-seeded and not cross-version/cross-language stable.

Each footgun affects *only* hash-table iteration order and runtime hash values.
Sorting the output and using content-defined keys makes every primitive a pure
function of its input multiset, immune to all three.

### Conformance-suite additions
`CONFORMANCE_SPEC.md` currently asserts *semantic* equivalence for graphs and
tolerates "minor formatting differences"; the bit-for-bit guarantee for relational
index sets is **new spec to add**. Tests should feed identical mesh inputs to all
three implementations and assert **byte-identical serialized index sets and
identical dense-ID arrays**, including adversarial inputs: duplicate edges, reversed
edge orientation, and permuted input order (to prove order-independence).

## 6. Reference API shape (canonical-sort convention)

Identical semantics in all three languages; only syntax differs. Hash structures
are used purely to bucket/dedup; the emitted order always comes from a sort on the
full tuple.

```
skolem_edge(a, b)      = (a <= b) ? (a, b) : (b, a)        # canonical tuple, no hash
distinct(tuples)       = sort(unique(tuples))               # total order, dedup
rank(distinct_sorted)  = { t -> i (+1 if 1-based) for i, t in enumerate(sorted) }
join(left, right, key) = bucket right by key; probe; emit sorted by canonical key
group_agg(rows, ⊕)     = bucket by group key; ⊕ within bucket; emit sorted by key
```

Concrete sketches per language (`indexmap::IndexSet` + `sort_unstable` in Rust;
`Set`+`sort!` in Julia; `np.unique(axis=0)`/`np.lexsort`/`searchsorted` in Python)
were produced during the research pass and should live next to each binding's
`relational` module.

## 7. Net recommendation

- **Architecture:** parallel native implementations, mirroring the toolkit. No
  shared engine, no heavy DataFrame/SQL dependency.
- **Libraries (all already depended-on):** Julia stdlib `Dict`/`Set`/`sort`
  (+`OrderedCollections`); Rust `indexmap`; Python NumPy.
- **The real work:** write §5 into `CONFORMANCE_SPEC.md` and back it with
  adversarial cross-binding tests. The five primitives are small and stable; the
  determinism spec is the durable artifact.
- **Effort:** ~100 lines per binding plus the spec; in Rust and Julia much of it is
  consolidating patterns (`performance.rs`/`canonicalize.rs`, `graph.jl`) that
  already exist; Python is greenfield but NumPy covers it directly.

## 8. Case study: conservative regridding (`ConservativeRegridding.jl` / `EarthSciData.jl`)

First-order conservative (area-weighted) regridding is a clean fourth specialization
of the IR — the same formal object as the FVM stencil and the ESI emissions
contraction, with the index space being *cell-overlap pairs between two grids*. The
mapping below is **verified against the actual `JuliaGeo/ConservativeRegridding.jl`
and `EarthSciML/EarthSciData.jl` source**, not just the abstract formula.

**The operation.** `F_target[j] = (1/A_j)·Σ_i A_ij·F_source[i]`, with
`A_ij = area(src_i ∩ tgt_j)` and `A_j = Σ_i A_ij`. In `ConservativeRegridding.jl`
this is a `Regridder` built once — a `SparseArrays` matrix of **raw overlap areas
`A_ij`** with row-sums stored as `dst_areas` — and applied by
`mul!(dst, intersections, src); dst ./= dst_areas`: a sparse matrix-vector product
plus a normalization divide. Build-once / apply-many is explicit in the API.

**FAQ decomposition (verified piece by piece):**

| Piece | `ConservativeRegridding.jl` reality | IR mapping | Partition |
|---|---|---|---|
| overlap pairs `{(i,j):A_ij>0}` | STR-tree dual depth-first search (`dual_depth_first_search`) — not O(n·m) | data-derived index set via a **spatial (theta) join**, executed by a spatial-index physical operator | static |
| `A_ij` | `GeometryOps.intersection` (Foster–Hormann planar / Sutherland–Hodgman spherical) + `GeometryOps.area`, `Manifold`-selectable | the **`intersect_polygon` kernel leaf** + **`polygon_area` FAQ** (RFC §8.1) — mirroring the `intersection`/`area` split | static |
| `A_j = Σ_i A_ij` | sparse row-sums → `dst_areas` | **group-by-`j` `sum_product` FAQ** | static |
| apply `Σ_i A_ij·F_src[i]` | `LinearAlgebra.mul!(dst, intersections, src)` | **`sum_product` FAQ** over the overlap set (sparse mat-vec) | dynamic (if `F_src` is time-varying) |
| `/A_j` | `dst ./= dst_areas` (or folded into the matrix at build, `normalize=true`) | elementwise; foldable to build time | static fold or dynamic |

Four of the five pieces are native FAQs; the fifth — `A_ij` — is the kernel-factor
leaf, which is exactly why RFC §8.1 makes `intersection_area` a required op. The
build-once/apply-many split the package already enforces *is* the static/dynamic
partition (RFC §6.1).

**The two boundaries, confirmed by the source:**

1. **Spatial join ≠ equi-join.** The overlap set is built with an **STR-tree dual
   DFS** — a spatial-index-accelerated theta-join — confirming the IR's equi-only
   `join.on` does not cover it. The declarative form is either the
   bin-Skolem-equi-join idiom or a first-class spatial-join kind backed by an
   STR/R-tree physical operator (a planner concern). The IR expresses *that* there
   is an overlap join; the spatial index is the operator that executes it.
2. **The clip is an opaque kernel; the area is not.** `GeometryOps.jl` already
   separates `intersection` (the polygon clip) from `area`. The IR mirrors this:
   `intersect_polygon` is the kernel leaf (iterative clipping, robustness-critical,
   not expressible as a semiring aggregate), while `polygon_area` is an ordinary
   `sum_product` FAQ over the clipped ring (shoelace / Gauss–Green / spherical excess)
   — which makes the weight differentiable w.r.t. geometry in-formalism (RFC §8.1).
   Only the clip's conformance is tolerance-based (floating-point), and it carries the
   planar-vs-spherical `Manifold` flag. Implementation options for the clipping leaf
   across the three bindings are the subject of the companion report
   [*Implementing the `intersect_polygon` kernel across ESS bindings*](intersection-area-kernel-implementation.md).

**`EarthSciData.jl` (verified):**
- **Non-staggered grids** use `ConservativeRegridding.jl` directly — the regridder is
  built once at `DataSetInterpolator` init and reused every step → the decomposition
  above.
- **Temporal interpolation** is a multilinear-in-time blend of cached slices
  (`DataSetInterpolator`/`TemporalCache`, temporal-first then spatial) → another
  `sum_product` FAQ (weighted sum of two time-slice factors). Its NetCDF loading /
  slice caching is data-source plumbing, **outside** the formalism (analogous to ESI
  table loading).
- **Staggered grids** use a separate `InterpolatingRegridder` (B-spline
  interpolation), *not* conservative regridding — a different operator. It is also a
  `sum_product` interpolation FAQ, but over a B-spline stencil with its own
  interpolation-weight factor rather than `intersection_area`. A second EarthSciData
  operator with its own kernel.

**Net:** conservative regridding's numeric core (overlap-area `sum_product` apply +
normalization group-by + temporal-interp blend) is native to the IR and inherits the
differentiability / XLA-tracing properties discussed for any `sum_product` FAQ (the
apply is a sparse mat-vec / `segment_sum` over precomputed indices). Its two non-FAQ
pieces — the spatial-index join operator and the `intersection_area` geometry kernel
— are exactly what the IR *coordinates* rather than subsumes. ConservativeRegridding's
`Regridder` is, in effect, a hand-built materialization of this static partition;
emitting it as ESS IR would unify it with ESM/ESD/ESI under one evaluator.

## 9. References

- ESS repo (inspected `main`): `packages/EarthSciSerialization.jl/{Project.toml,src/tree_walk.jl,src/graph.jl}`,
  `packages/earthsci-toolkit-rs/{Cargo.toml,src/performance.rs,src/canonicalize.rs}`,
  `packages/earthsci_toolkit/{pyproject.toml,src/earthsci_toolkit/numpy_interpreter.py,canonicalize.py}`,
  `CONFORMANCE_SPEC.md` — <https://github.com/EarthSciML/EarthSciSerialization>.
- Determinism: Rust `HashMap` SipHash randomization (`std::collections::HashMap`
  docs); PEP 456 (Python SipHash / hash randomization); `indexmap` hasher-independent
  order (docs.rs/indexmap); Julia `sort` stability (≥1.9).
- Cross-language hashing: xxHash / XXH3 (`xxhash-rust`, `xxhash` PyPI, `XXhash.jl`).
- Engines surveyed: DuckDB 1.5.2 (`DuckDB.jl`, `duckdb`, `duckdb-rs`); Polars
  (`polars` crate + PyPI; `Polars.jl` unmaintained); Apache Arrow / Acero
  (`pyarrow` only; `Arrow.jl` format-only; Rust DataFusion).
- Parent RFC: [`semiring-faq-unified-ir.md`](semiring-faq-unified-ir.md) §6.1, §9.
</content>
</invoke>
