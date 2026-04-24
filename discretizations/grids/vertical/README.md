# Vertical grid family — canonical fixtures

Declarative `.esm` grid configs for the `vertical` family, at standard
resolutions covering all coordinate kinds supported by
`docs/GRIDS_API.md` §2.4 / §3 (see also the Python binding at
`python/src/earthsci_toolkit/grids/vertical.py`).

Per the 2026-04-20 mayor correction on bead `dsc-0he`, each `.esm` file is a
small **declarative config** — family + coordinate kind + interface levels
(+ hybrid A/B coefficients when applicable) + provenance — not a serialized
geometry blob. Centres and widths are derived on demand by each binding's
accessor runtime.

## Layout

```
discretizations/grids/
├── vertical.schema.json           # family schema (peer of cartesian/lat_lon)
└── vertical/
    ├── README.md                  # this file
    ├── regenerate_fixtures.py     # generator (reference-binding = Python)
    ├── sigma_uniform_n16.esm      # uniform sigma, nz = 16
    ├── sigma_uniform_n64.esm      # uniform sigma, nz = 64
    ├── z_troposphere_l32.esm      # geometric altitude, 32 layers to ~34 km
    ├── eta_hybrid_l12.esm         # CAM-style hybrid A/B, 12 layers
    └── theta_isentropic_l10.esm   # isentropic theta, 280 K → 380 K
```

## Fixture provenance

| Fixture | Coordinate | `nz` | Generator | Config summary |
|---|---|---|---|---|
| `sigma_uniform_n16`   | `sigma` | 16 | Python binding, `vertical(coordinate="sigma", nz=16)` | Uniform sigma; interface 0 = surface (σ=1), interface nz = top (σ=0). Small resolution for fast tests. |
| `sigma_uniform_n64`   | `sigma` | 64 | Python binding, `vertical(coordinate="sigma", nz=64)` | Uniform sigma, medium resolution. |
| `z_troposphere_l32`   | `z`     | 32 | Python binding, explicit `levels` | Increasing-thickness geometric altitude from 0 m to 34 km. Typical research-model stretched grid for the troposphere + lower stratosphere. |
| `eta_hybrid_l12`      | `eta`   | 12 | Python binding, `ak` / `bk` table | Hybrid sigma-pressure, CAM-style A/B shape. `p0 = 1.0e5 Pa`. Synthesized σ = `ak/p0 + bk` is strictly decreasing from 1.0 at the surface to 0.0 at the model top. |
| `theta_isentropic_l10`| `theta` | 10 | Python binding, explicit `levels` | Isentropic potential-temperature surfaces from 280 K to 380 K in 10 K steps. |

The Python binding (`python/src/earthsci_toolkit/grids/vertical.py`) is the
**reference binding** for these fixtures today because it is the first
implementation to land (per `dsc-g4u`). As the Julia (`dsc-79n`), Rust
(`dsc-31p`), and TypeScript (`dsc-0k6`) accessor runtimes land and the
vertical conformance harness (`dsc-bbj`) comes online, each binding will be
required to emit byte-identical declarative `.esm` content when called with
the same options.

## Reference binding rationale

Per `docs/GRIDS_API.md` §4.3, Julia is the conformance reference binding for
ULP ties. For these declarative fixtures there is no transcendental math —
sigma levels are pure rationals and geometric / theta / eta profiles are
stored verbatim — so every binding should reproduce the same bytes
regardless of libm. Python is used as the *generator of record* here purely
because it landed first; the conformance bead will replace the generator
with the reference-binding (Julia) output on landing.

## Regenerating

```bash
# From the repo root
PYTHONPATH=python/src python3 discretizations/grids/vertical/regenerate_fixtures.py
```

The script is idempotent and deterministic. Running it against a clean
working tree should produce no diffs. Commit only intentional changes; the
Julia inline tests in `test/test_vertical_fixtures.jl` will fail CI if a
fixture drifts out of schema (missing `provenance`, wrong monotonicity,
inconsistent hybrid coefficients, …).

## CI coverage

The Julia TestItemRunner suite `test/test_vertical_fixtures.jl` executes on
every PR via `Pkg.test()`. It verifies:

1. `vertical.schema.json` is present and well-formed.
2. At least two `.esm` fixtures exist in this directory.
3. Each fixture parses and has the `docs/GRIDS_API.md` §7 common-minimum
   shape (`family`, `topology`, `ndim`, `dtype`, `ghosts`, `n_cells`,
   `n_vertices`, `n_edges`, `schema_version`, `provenance`).
4. `levels` monotonicity matches the coordinate kind (strictly decreasing
   for `sigma` / `eta` / `hybrid_sigma_theta`, strictly increasing
   otherwise).
5. Hybrid fixtures (`eta`) carry `ak`, `bk`, `p0`, with `ak.size == bk.size
   == nz + 1` and synthesized σ strictly decreasing.
6. The canonical corpus covers `sigma` (≥ 2 resolutions), `z`, `eta`,
   `theta`.

This is the "Layer A" byte-level schema conformance called out by the ESD CI
walker (`dsc-usq`). The cross-binding byte-equality check (comparing these
declarative bytes against Julia / Rust / TypeScript `to_esm` output) will be
wired up in the vertical conformance bead (`dsc-bbj`).
