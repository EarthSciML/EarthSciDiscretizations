# Projected (Lambert Conformal Conic) grid construction

Declarative `.esm` documents for **projected native grids** in the GDD Grid
`crs` family (RFC pure-io-data-loaders §4.2 + §5.3, bead `esd-47z.5`). A
projected grid is topologically `cartesian` — a uniform lattice in projected
metres — carrying a [`crs` descriptor](../cartesian.schema.json) (projection,
datum, radius, parameters) that names the projection a downstream rule consumes.
The `crs` descriptor is **orthogonal to the topological family** (ess-v9a.2): a
lat-lon grid and a Lambert-Conformal grid both have `family: cartesian` and
differ only in `crs`.

**Constructing** the grid's geographic geometry means recovering the
`(lon, lat)` of every cell **corner** by applying the spherical Lambert
Conformal Conic **inverse** to the regular projected `(x, y)` corner lattice —
closed-form, no iteration. The inverse is the inverse half of the reprojection
family member [`reprojection/lambert_conformal.esm`](../../../reprojection/lambert_conformal.esm)
(`esd-47z.2`); this construction reuses its cone constants and inverse formulas,
applied to lattice corners instead of round-tripped forward points.

## Files

| File | Role |
|------|------|
| `wrf_lcc.esm` | WRF native LCC projected Grid descriptor (`crs` + lattice). |
| `nei2016_lcc.esm` | NEI2016 native LCC projected Grid descriptor. |
| `rules/lambert_conformal_construction.esm` | The construction rule: spherical-LCC inverse over a projected corner `(x, y)` → cell-corner `(lon, lat)`. |
| `lambert_conformal_construction_conformance_test.jl` | Julia evaluator conformance (round-trip + native geometry). |

## Inverse (Snyder 1987, eqs. 15-5..15-9 = PROJ `+proj=lcc +R=…`)

With cone constants fixed by the `crs` parameters
`n = log(cos φ₁/cos φ₂)/log(tan(π/4+φ₂/2)/tan(π/4+φ₁/2))`,
`F = cos φ₁·tan(π/4+φ₁/2)^n/n`, `ρ₀ = R·F/tan(π/4+φ₀/2)^n`:

```
ρ   = sign(n)·sqrt(x² + (ρ₀−y)²)
θ   = atan2(x, ρ₀−y)
lon = lon_0 + (θ/n)·180/π
lat = (2·atan((R·F/ρ)^(1/n)) − π/2)·180/π
```

Like the rest of the `reprojection/` family it is a **scalar point-wise kernel**
(`proj_x`/`proj_y` scalar `parameter`s; `geo_lon`/`geo_lat` scalar `observed`s),
authored over the existing operator set (`sin`/`cos`/`tan`/`atan`/`atan2`/`log`/
`^`/`sqrt`/`sign` + the built-in `pi` op) with **no** new rule-schema, scaffold,
validation, resolver, or per-binding code. The recovered corner coordinates are
surfaced through `build_evaluator` via the family's zero-IC state idiom
(`d(corner_lon)/dt = geo_lon` from a zero IC ⇒ `du = geo_lon`). The conformance
harness sweeps the kernel over the `(Nx+1)×(Ny+1)` corner lattice each Grid
descriptor pins (a binding applies the point map over the corner array at the
runtime layer).

## Parameter sets (one construction, two real native grids)

| Grid | `lat_1` | `lat_2` | `lat_0` | `lon_0` | `R` (m) | cone `n` |
|------|---------|---------|---------|---------|---------|----------|
| WRF | 30 | 60 | 38.999996 | −97 | 6 370 000 | ≈ 0.7156 |
| NEI2016 | 33 | 45 | 40 | −97 | 6 370 997 | ≈ 0.6305 |

Both grids pin the same representative CONUS lattice (`x0=−2e6`, `dx=1e6`,
`Nx=4`; `y0=−1.5e6`, `dy=1e6`, `Ny=3` ⇒ 5×4 corners spanning lon∈[−126°,−68°],
lat∈[23°,53°], well away from the singular far pole). The cone constant depends
only on the standard parallels, so the same lattice recovers different geography
under the two `crs` parameter sets — the declarative-or-fail proof the
construction binds a genuinely **parameterized** grid.

## Acceptance (`esd-47z.5`)

* Both projected Grid descriptors **validate and load** with the `crs` descriptor
  preserved (the GDD Grid crs round-trip).
* The constructed cell-corner `(lon, lat)` match an **independent proj4/Snyder
  reference** (computed in NumPy outside the `.esm`) for both grids.
* **forward∘inverse = identity**: forward-projecting the constructed corners
  through `reprojection/lambert_conformal.esm` restores the projected `(x, y)`
  lattice (the two reprojection-family members compose to the identity).
* The `crs` parameters are **load-bearing**: WRF ≠ NEI2016 off the shared
  central meridian.
