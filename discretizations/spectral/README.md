# Spectral Rules

Rule files for spectral and pseudo-spectral discretizations (Fourier,
Chebyshev, spherical-harmonic transforms, etc.).

The convention for file naming and organization within this directory may
evolve as content lands. Current expectation: one JSON file per named scheme,
e.g. `fourier_1d.json`, `spherical_harmonic_t63.json`.

## Status: no catalog rules yet — see the feasibility verdict

The declarative feasibility of this family was investigated under bead
esd-6g4.4 (DECLARATIVE-OR-FAIL). Verdict: **`SPECTRAL_FEASIBILITY.md`**.

Summary: spectral operators are global, expressible only as a dense
`reduce`-`arrayop` contraction. That mechanism works, and the Fourier and
Chebyshev differentiation matrices are fully AST-expressible (proven to machine
precision — see `feasibility_probe/`). But a **general** rule is blocked by an
ESS gap (a reduction-range bound cannot track the grid dimension size), and
spherical harmonics are fundamentally infeasible (associated Legendre functions
are recurrence-defined; the AST has no recurrence). No imperative solver was
added. Read the verdict before authoring spectral rules.
