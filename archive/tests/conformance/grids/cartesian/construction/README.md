# Cartesian construction value-FAQ conformance

Cross-binding byte-identity harness for **cartesian structured-grid construction**
expressed as a declarative elementwise FAQ (bead `esd-3we.1` / S1). This is the
**structured-grid FAQ template** that the lat-lon / vertical / Arakawa steps
(S2–S4) reuse — the structured-grid counterpart to the DUO
`grids/duo/topology/` value-invention harness. It pins, for the shared cartesian
conformance fixtures, the construction arrays the landed **M1** elementwise FAQ
produces from each fixture's grid parameters:

| array              | what it is                                          | M1 elementwise op |
|--------------------|-----------------------------------------------------|-------------------|
| `edges`            | per-axis cell faces (`lo + (i-1)*dx`, or supplied)  | affine map / data |
| `centers`          | per-axis cell centers (`(e[i]+e[i+1])/2`)           | elementwise midpoint |
| `widths`           | per-axis cell widths (`e[i+1]-e[i]`)                | elementwise difference |
| `cell_volume`      | per-cell measure (∏ of the axis widths)             | elementwise product |
| `metric_jacobian`  | identity-metric Jacobian (= `cell_volume`)          | identity metric |
| `neighbor_minus/plus` | offset ∓1 axis neighbor, `-1` = no neighbor      | index arithmetic (0-sentinel) |
| `boundary_lower/upper` | first / last cell of the axis (`0`/`1`)          | boundary mask |

## Artifacts

- **`../../../../../discretizations/grids/cartesian/rules/cartesian_construction.esm`**
  — the declarative M1 elementwise FAQ document (RFC `semiring-faq-unified-ir`
  §5.1/§5.2). Declares the affine coordinates, the identity metric
  (`metric_g`/`metric_ginv` = `delta_ij`, `jacobian` = ∏ widths, `dgij_dxk` = 0,
  `coord_jacobian` = I, `coord_jacobian_second` = 0), the 0-sentinel neighbor
  maps, the boundary masks and `cell_volume` over two interval index sets
  (`edge_nodes`, `cells`).
- **`src/cartesian_faq.jl` → `cartesian_construction_faq`** — the thin evaluation
  bridge that routes the affine map + every coordinate/volume arithmetic step
  through `EarthSciSerialization`'s AST evaluator (`eval_coeff`). ESD hosts no
  shadow evaluator (AGENTS.md single-pathway; GRIDS_API §4.3); the structural
  arrays (cartesian-product flattening, neighbor linearization, masks, identity
  constants) are pure index logic in the bridge.
- **`golden.json`** — this directory. The binding-neutral serialization of the
  construction arrays for every fixture in `../fixtures.json`.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge, reading
  the SHARED `../fixtures.json`.

## Serialization & indexing (binding-neutral)

Each construction array family is stored as **one compact JSON string** of its
`[axis][…]` nesting — the same dense byte form as
`../../duo/topology/golden.json`. Julia is the reference binding (1-based ids,
`0` = no neighbor):

- **Coordinate / metric floats** (`edges`, `centers`, `widths`, `cell_volume`,
  `metric_jacobian`) are full-precision `Float64`. The Julia reference test
  compares the strings **exactly** (the construction is pure multiplication and
  linear interpolation — deterministic IEEE-754, deterministic Ryu repr — so the
  bytes are stable). The **cross-binding** contract compares these at the family
  relative tolerance (`tolerance.relative`, `1e-14`), since another binding's
  float formatting may differ while the value does not.
- **Neighbor maps** (`neighbor_minus`, `neighbor_plus`) are **0-based** with
  boundary sentinel **`-1`** (Julia's 1-based linear ids and `0` = no-neighbor
  converted at the harness boundary, identical to `../regenerate_golden.jl`);
  **boundary masks** are `0`/`1`. Integer arrays match exactly in every binding.

## Why it is byte-identical across Julia / Rust / Python

The grid parameters (`extent`, `n`, supplied `edges`) are the same in every
binding, and the construction arithmetic — the affine map and the
midpoint/difference/product derivations — is evaluated by `EarthSciSerialization`'s
AST evaluator, whose result is a pure function of the IEEE-754 operations in a
fixed order (no transcendentals, never hash or insertion order), guaranteed
identical across the three bindings by the ESS determinism contract
(`CONFORMANCE_SPEC.md` §5.5). So each binding's `cartesian_construction_faq` output
reproduces these values, and the integer neighbor / boundary arrays match exactly.

## What the test asserts (`test/test_cartesian_construction_faq.jl`)

For each fixture (`small_1d`, `small_2d`, `realistic_3d`, `nonuniform_2d`):
1. **Match imperative cartesian.jl to ULP** — every FAQ array (coords, identity
   metric, neighbor maps, boundary masks) equals the imperative
   `src/grids/cartesian.jl` trait array bit-for-bit.
2. **Internal consistency** — edges ascend, widths are positive, each center lies
   strictly inside its cell, the metric is the Kronecker delta with vanishing
   derivatives, and the neighbor maps obey the 0-sentinel and ∓ symmetry.
3. **Byte-identity** — the binding-neutral serialization equals `golden.json`.

The `.esm` document's schema validity is also asserted.

## Regenerate

```bash
julia --project=. tests/conformance/grids/cartesian/construction/regenerate_golden.jl
```

Regenerate only when the construction contract legitimately changes; an unexpected
diff means the FAQ output drifted from the imperative builder.
