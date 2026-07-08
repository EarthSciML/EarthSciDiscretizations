# Arakawa construction value-FAQ conformance

Cross-binding byte-identity harness for **arakawa staggered-grid construction**
expressed as a declarative elementwise FAQ (bead `esd-3we.4` / S4). The staggered
sibling of the cartesian `../../cartesian/construction/` structured-grid FAQ template.
It pins, for the shared arakawa conformance fixtures, the construction the landed
**M1** elementwise FAQ produces from each fixture's base-grid parameters:

| output                | what it is                                              | M1 op |
|-----------------------|---------------------------------------------------------|-------|
| `center` map          | per-axis cell centres (`lo + (i-0.5)*d`)                | affine map (½-cell stagger) |
| `node` map            | per-axis cell faces / corners (`lo + (i-1)*d`)          | affine map |
| `cell_center/u_edge/v_edge/corner` | the four staggered locations          | cartesian product of the two maps |
| `cell_volume`         | per-cell area (`dx*dy`)                                  | elementwise product |
| `neighbor_*`          | row-major ∓1 axis neighbor, `-1` = no neighbor          | index arithmetic (0-sentinel) |
| `boundary_*`          | first / last cell of the axis (`0`/`1`)                 | boundary mask |
| `h/u/v_location`      | static A–E variable-location table                      | closed-form lookup (`ifelse`) |

The arakawa stagger over a cartesian base is **separable**: each of the four
locations is the cartesian product of two staggered 1-D affine maps per axis — a
`center` map (cell centres, half-cell offset) and a `node` map (cell faces /
corners). The `node` map is exactly the cartesian-template `edges` map, so the ESS
evaluation is bit-identical to the imperative builder; the `center` map adds only the
`i - 0.5` half-cell offset.

## Artifacts

- **`../../../../../discretizations/grids/arakawa/rules/arakawa_construction.esm`** —
  the declarative M1 elementwise FAQ document (RFC `semiring-faq-unified-ir`
  §5.1/§5.2). Declares the staggered `center`/`node` axis maps, the `dx*dy`
  cell-volume product, the 0-sentinel neighbor maps, the boundary masks, and the
  static A–E variable-location table (`h/u/v_location` over a `staggers` interval).
- **`src/arakawa_faq.jl` → `arakawa_construction_faq`** — the thin evaluation bridge
  that routes `dx`/`dy`, the two staggered axis maps, and the cell-volume product
  through `EarthSciSerialization`'s AST evaluator (`eval_coeff`). ESD hosts no shadow
  evaluator (AGENTS.md single-pathway; GRIDS_API §4.3); the structural arrays (the
  cartesian-product flattening into the four staggered locations, the row-major
  neighbor linearization, the masks, and the static location table) are pure index
  logic in the bridge.
- **`golden.json`** — this directory. The binding-neutral serialization of the
  construction for every fixture in `../fixtures.json`.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge, reading the
  SHARED `../fixtures.json`.

## Serialization & indexing (binding-neutral)

The shared arakawa fixtures are large (the `realistic` base is 256×256), so — like
the accessor golden `../golden/*.json` — the coordinate / neighbor / boundary arrays
are pinned at the fixtures' `query_points` rather than in full. (The Julia test still
asserts the **full** arrays match the imperative builder to ULP in memory; only the
byte golden is sampled.) The static A–E table is small and pinned in full per stagger.

- **Coordinate / metric floats** (`dx`, `dy`, `cell_volume`, the per-location `xy`)
  are full-precision `Float64`. The Julia reference test compares the compact
  coordinate strings **exactly** (pure multiply/add — deterministic IEEE-754,
  deterministic Ryu repr); the **cross-binding** contract compares these at the
  family relative tolerance (`tolerance.relative`, `1e-14`).
- **Neighbor maps** are **0-based** with boundary sentinel **`-1`** (Julia's 1-based
  row-major linear ids and `0` = no-neighbor converted at the harness boundary);
  **boundary masks** are `0`/`1`; **location codes** are the names
  `cell_center`/`u_edge`/`v_edge`/`corner`. Integer arrays match exactly in every
  binding.

## Why it is byte-identical across Julia / Rust / Python

The base-grid parameters (`extent`, `nx`, `ny`) are the same in every binding, and
the construction arithmetic — the two affine maps and the product — is evaluated by
`EarthSciSerialization`'s AST evaluator, whose result is a pure function of the
IEEE-754 operations in a fixed order (no transcendentals, never hash or insertion
order), guaranteed identical across the three bindings by the ESS determinism
contract (`CONFORMANCE_SPEC.md` §5.5). So each binding's `arakawa_construction_faq`
output reproduces these values, and the integer neighbor / boundary / location-table
arrays match exactly.

## What the test asserts (`test/test_arakawa_construction_faq.jl`)

For each fixture (`small`, `realistic`) × each stagger (`A`–`E`):
1. **Match imperative arakawa.jl to ULP** — every FAQ array (the four staggered
   location coordinates, the Tier-C cell-centred widths/volume, the row-major neighbor
   maps and boundary masks, the metric `dx`/`dy`/`area`) equals the imperative
   `src/grids/arakawa.jl` accessor bit-for-bit, and the static A–E
   variable-location / shape table equals `arakawa_variable_locations` /
   `variable_shape` / `arakawa_shape`.
2. **Internal consistency** — nodes ascend, each centre lies strictly inside its
   cell, widths are positive, `cell_volume = dx*dy`, and the row-major neighbor maps
   obey the 0-sentinel and ∓ symmetry.
3. **Byte-identity** — the binding-neutral serialization equals `golden.json`.

The `.esm` document's schema validity is also asserted.

## Regenerate

```bash
julia --project=. tests/conformance/grids/arakawa/construction/regenerate_golden.jl
```

Regenerate only when the construction contract legitimately changes; an unexpected
diff means the FAQ output drifted from the imperative builder.
