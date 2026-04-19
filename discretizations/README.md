# Discretization Rule Files

This directory holds authoritative discretization rule JSON files. Each file
describes how a continuous PDE operator maps onto a discrete stencil, grid
staggering, and boundary-condition handling.

Rule files are validated against the EarthSciSerialization discretization
schema (§7). Validation and application are performed by the ESS rule engine
once it lands; until then, these files are parsed as opaque JSON.

## Layout

- `finite_difference/` — finite-difference stencils on structured grids
- `finite_volume/` — finite-volume reconstructions and flux forms
- `spectral/` — spectral / pseudo-spectral methods

The per-family split is a starting convention and may evolve as content lands.
See `../docs/REPO_LAYOUT.md` for the governing convention.

## Adding a rule

1. Drop a `.json` file into the appropriate family subdirectory
2. Ensure it validates against the ESS §7 schema
3. Add a rule-application test under `../test/` (runs against the rule engine)
