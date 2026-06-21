# Vertical construction value-FAQ conformance

Cross-binding byte-identity harness for **vertical (1-D column) structured-grid
construction** expressed as a declarative elementwise FAQ (bead `esd-3we.2` / S2).
It **reuses the structured-grid FAQ template** established by cartesian
(`../../cartesian/construction/`, `esd-3we.1` / S1) — the 1-D column counterpart.
It pins, for the shared vertical conformance fixtures, the construction arrays the
landed **M1** elementwise FAQ produces from each fixture's column parameters:

| array                 | what it is                                              | M1 elementwise op |
|-----------------------|--------------------------------------------------------|-------------------|
| `levels`              | interface (half-level) coords, synthesized per kind    | affine / multiply-add / data |
| `centers`             | cell centers (`(levels[k]+levels[k+1])/2`)             | midpoint reduction |
| `widths`              | layer thicknesses (`abs(levels[k+1]-levels[k])`)       | absolute difference |
| `cell_volume`         | per-layer measure (= `widths` in 1-D)                  | identity |
| `metric_dz` / `metric_z` | thickness / cell-center metric                      | identity read |
| `metric_sigma`        | sigma at center (sigma-like kinds)                     | identity read |
| `metric_pressure`     | `((ak+bk·p0)+(ak'+bk'·p0))/2` layer-average            | layer-average |
| `metric_ak` / `metric_bk` | coefficient layer-averages (`(c[k]+c[k+1])/2`)     | layer-average |
| `neighbor_minus/plus` | offset ∓1 layer neighbor, `-1` = no neighbor           | index arithmetic (0-sentinel) |
| `boundary_lower/upper`| bottom / top layer of the column (`0`/`1`)             | boundary mask |
| `ak` / `bk` / `p0`    | `pressure_coefficients` (ak/bk present for hybrid)     | data |

## Level synthesis per coordinate kind

- **uniform sigma** (`sigma`, `hybrid_sigma_theta` without supplied levels):
  `levels[k] = 1 - (k-1)/nz` — the affine-like analogue of cartesian's
  `edges[k] = lo + (k-1)·dx`.
- **eta** hybrid sigma-pressure: `levels[k] = ak[k]/p0 + bk[k]`, where `ak`/`bk`
  are DATA inputs — one divide-add per interface.
- **explicit** `z` / `theta` / `z_star` (and supplied-`levels` sigma): the
  interface coordinates are supplied verbatim as DATA (the supplied-array branch,
  the analogue of cartesian's non-uniform `edges`).

## Artifacts

- **`../../../../../discretizations/grids/vertical/rules/vertical_construction.esm`**
  — the declarative M1 elementwise FAQ document (RFC `semiring-faq-unified-ir`
  §5.1/§5.2). Two models over the `interfaces` (nz+1) / `layers` (nz) interval
  index sets: `VerticalSigmaConstruction` (pinned to `sigma_uniform_n16`: uniform
  level synthesis, centers/widths/volume, dz/z/sigma metrics, neighbor maps,
  boundary masks) and `VerticalHybridConstruction` (pinned to `eta_hybrid_l12`: the
  `ak/p0+bk` synthesis and the pressure/ak/bk layer-average metrics).
- **`src/vertical_faq.jl` → `vertical_construction_faq`** — the thin evaluation
  bridge that routes every level/center/width/metric arithmetic step through
  `EarthSciSerialization`'s AST evaluator (`eval_coeff`). ESD hosts no shadow
  evaluator (AGENTS.md single-pathway; GRIDS_API §4.3); the structural integer
  arrays (neighbor maps, masks) are pure index logic in the bridge.
- **`golden.json`** — this directory. The binding-neutral serialization of the
  construction arrays for every fixture in `../fixtures.json`.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge, reading
  the SHARED `../fixtures.json`.

## Serialization & indexing (binding-neutral)

Each construction array is stored as **one compact JSON string** (no spaces) — the
same dense byte form as the cartesian construction golden. Julia is the reference
binding (1-based ids, `0` = no neighbor):

- **Level / coordinate / metric floats** (`levels`, `centers`, `widths`,
  `cell_volume`, `metric_*`, `ak`, `bk`) are full-precision `Float64`. The Julia
  reference test compares the strings **exactly**; the vertical family has no
  transcendental math (pure rationals + a single multiply-add), so the family
  tolerance is `0.0` (strict byte equality) — the **cross-binding** contract also
  compares these exactly.
- **Neighbor maps** (`neighbor_minus`, `neighbor_plus`) are **0-based** with
  boundary sentinel **`-1`** (Julia's 1-based ids and `0` = no-neighbor converted at
  the harness boundary, identical to the cartesian construction golden);
  **boundary masks** are `0`/`1`. Integer arrays match exactly in every binding.

## Why it is byte-identical across Julia / Rust / Python

The column parameters (`coordinate`, `nz`, supplied `levels`, hybrid `ak`/`bk`/`p0`)
are the same in every binding, and the construction arithmetic — the level
synthesis, the midpoint/difference, and the layer-average metrics — is evaluated by
`EarthSciSerialization`'s AST evaluator, whose result is a pure function of the
IEEE-754 operations in a fixed order (no transcendentals, never hash or insertion
order), guaranteed identical across the three bindings by the ESS determinism
contract (`CONFORMANCE_SPEC.md` §5.5). So each binding's `vertical_construction_faq`
output reproduces these values, and the integer neighbor / boundary arrays match
exactly.

## What the test asserts (`test/test_vertical_construction_faq.jl`)

For each fixture (`sigma_uniform_n16`, `sigma_uniform_n64`, `z_troposphere_l32`,
`eta_hybrid_l12`, `theta_isentropic_l10`):
1. **Match imperative vertical.jl to ULP** — every FAQ array (levels, centers,
   widths, volume, the available metrics via `metric_eval`, neighbor maps, boundary
   masks, pressure_coefficients) equals the imperative `src/grids/vertical.jl` trait
   array / accessor bit-for-bit.
2. **Internal consistency** — `length(levels) == nz+1`, strictly monotone levels,
   positive widths, each center strictly between its interfaces, the neighbor maps
   obey the 0-sentinel and ∓ symmetry, and `half_levels == levels` /
   `layer_thickness == widths`.
3. **Byte-identity** — the binding-neutral serialization equals `golden.json`.

The `.esm` document's schema validity is also asserted.

## Regenerate

```bash
julia --project=. tests/conformance/grids/vertical/construction/regenerate_golden.jl
```

Regenerate only when the construction contract legitimately changes; an unexpected
diff means the FAQ output drifted from the imperative builder.
