# Lat-lon construction value-FAQ conformance

Cross-binding byte-identity harness for **lat-lon structured-grid construction**
expressed as a declarative elementwise FAQ (bead `esd-3we.3` / S2). This is the
lat-lon instance of the **structured-grid FAQ template** established by cartesian
(`esd-3we.1` / S1) — the rectilinear-sphere counterpart to the DUO
`grids/duo/topology/` value-invention harness. It pins, for the shared lat-lon
conformance fixtures (regular + reduced-Gaussian), the construction arrays the landed
**M1** elementwise FAQ produces from each fixture's grid parameters:

| array                       | what it is                                                   | M1 elementwise op |
|-----------------------------|-------------------------------------------------------------|-------------------|
| `cell_centers_lon/lat`      | per-cell `(lon, lat)` (`lon_start + (i-½)·dlon`; row center) | affine map / data |
| `cell_widths_lon/lat`       | per-cell `dlon` and the lat-edge difference                 | affine / difference |
| `cell_volume`               | spherical-rectangle area `R²·dlon·(sin φ_n − sin φ_s)`       | trig measure |
| `lat_edges`, `lat_centers`  | 1-D latitude faces (data) and row centers (midpoint/data)   | data / midpoint |
| `metric_g_lonlon`           | `g_λλ = R² cos²φ`                                            | trig metric |
| `metric_g_latlat`           | `g_φφ = R²`                                                  | constant metric |
| `metric_ginv_lonlon`        | `g^λλ = 1/(R² cos²φ)` (Inf at the poles)                     | trig metric inverse |
| `metric_jacobian`           | `J = R² \|cos φ\|`                                           | trig Jacobian |
| `dg_lonlon_dlat`            | `∂g_λλ/∂φ = −2R² cos φ sin φ` (the only non-zero derivative) | trig derivative |
| `neighbor_lon_minus/plus`   | periodic-longitude wrap (no sentinel)                       | index arithmetic |
| `neighbor_lat_minus/plus`   | nearest-center adjacent row, `-1` = pole (`pole_policy=:none`) | rank/join (0-sentinel) |
| `boundary_lat_lower/upper`  | first / last latitude row (`0`/`1`)                         | boundary mask |

## Artifacts

- **`../../../../../discretizations/grids/latlon/rules/latlon_construction.esm`**
  — the declarative M1 elementwise FAQ document (RFC `semiring-faq-unified-ir`
  §5.1/§5.2). Declares the latitude coords (data + midpoint/difference), the per-row
  affine longitude map, the spherical-rectangle area, the closed-form curvilinear
  metric (`g_λλ`/`g_φφ`/`g_λφ`, `ginv`, `J`, `∂g_λλ/∂φ`, `coord_jacobian` = I,
  `coord_jacobian_second` = 0), the periodic / nearest-center neighbor maps and the
  latitude boundary masks over four interval index sets (`lat_nodes`, `rows`,
  `lon_cells`, `cells`).
- **`src/latlon_faq.jl` → `latlon_construction_faq`** — the thin evaluation bridge
  that routes the affine longitude map, the cell area, and every metric component
  (`sin`/`cos`/`abs`) through `EarthSciSerialization`'s AST evaluator (`eval_coeff`).
  ESD hosts no shadow evaluator (AGENTS.md single-pathway; GRIDS_API §4.3); the
  structural arrays (ragged-row flattening, periodic / nearest-center neighbor
  linearization, masks, identity coordinate Jacobian) are pure index logic in the
  bridge — the lat-lon analogue of `cartesian_faq.jl`'s neighbor linearization.
- **`golden.json`** — this directory. The binding-neutral serialization of the
  construction arrays for every fixture in `../fixtures.json`.
- **`regenerate_golden.jl`** — regenerates `golden.json` from the bridge, reading the
  SHARED `../fixtures.json`.

## Serialization & indexing (binding-neutral)

Each construction array is stored as **one compact JSON string** of its flat
ragged-row-major cell vector (or the 1-D latitude vector). Julia is the reference
binding for this golden (1-based ids, `0` = no neighbor):

- **Coordinate / metric floats** (`cell_centers_*`, `cell_widths_*`, `cell_volume`,
  `lat_edges`, `lat_centers`, `metric_*`, `dg_lonlon_dlat`) are full-precision
  `Float64`. The Julia reference test compares the strings **exactly**; the
  **cross-binding** contract compares them at the family relative tolerance
  (`tolerance.relative`, `1e-14`), since the latitudes feed `sin`/`cos` and another
  binding's libm may drift by sub-ULP amounts while another binding's float formatting
  may differ.
- **Neighbor maps** (`neighbor_lon_minus/plus`, `neighbor_lat_minus/plus`) are
  **0-based** with pole sentinel **`-1`** (Julia's 1-based flat ids and `0` = no
  neighbor converted at the harness boundary, identical to
  `../../cartesian/construction/regenerate_golden.jl`); **boundary masks** are `0`/`1`.
  Integer arrays match exactly in every binding. Longitude neighbors wrap periodically
  (never `-1`); latitude neighbors are `-1` only at the poles.

## Why it is byte-identical across Julia / Rust / Python

The grid parameters (`R`, `nlon`/`nlon_per_row`, `lat_edges`, `lon_start`) are the same
in every binding, and the construction arithmetic — the affine longitude map, the
spherical-rectangle area, and the trig metric — is evaluated by
`EarthSciSerialization`'s AST evaluator, whose result is a pure function of the
IEEE-754 operations (including `sin`/`cos`/`abs`) in a fixed order, guaranteed
identical across the three bindings by the ESS determinism contract
(`CONFORMANCE_SPEC.md` §5.5). So each binding's `latlon_construction_faq` output
reproduces these values, and the integer neighbor / boundary arrays match exactly.

## What the test asserts (`test/test_latlon_construction_faq.jl`)

For each fixture (`small`, `realistic`, `reduced`):
1. **Match imperative latlon.jl to ULP** — every FAQ array (lon/lat coords, widths,
   area, trig metric, neighbor maps, boundary masks) equals the imperative
   `src/grids/latlon.jl` trait array bit-for-bit, including the trig metric, the
   periodic longitude neighbors, and the reduced-Gaussian nearest-center remap.
2. **Internal consistency** — areas positive, the orthogonal-metric identities
   (`g_λφ = 0`, `g_λλ·g^λλ = 1`, `J = R²|cos φ|`), the periodic-longitude `±` symmetry,
   and the latitude pole sentinel.
3. **Byte-identity** — the binding-neutral serialization equals `golden.json`.

The `.esm` document's schema validity is also asserted.

## Regenerate

```bash
julia --project=. tests/conformance/grids/latlon/construction/regenerate_golden.jl
```

Regenerate only when the construction contract legitimately changes; an unexpected
diff means the FAQ output drifted from the imperative builder.
