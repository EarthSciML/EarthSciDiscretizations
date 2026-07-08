# Reprojection cross-binding conformance harness

Verifies the declarative reprojection models (`reprojection/<model>.esm`) are
green **through the ESS engine** and pins the binding-neutral fixtures + golden
that every binding (Julia, Python, Rust, TypeScript) consumes.

Scope follows RFC `pure-io-data-loaders` §8 and `CONFORMANCE_SPEC.md` §5.8
(bead `esd-47z.6`, the "ESD DRAIN" conformance bead for the
`reprojection/`+`regridding/` family). The per-model **Julia evaluator** tests
ship co-located with each `.esm`
(`reprojection/{longlat,lambert_conformal}_conformance_test.jl`, beads
`esd-47z.1`/`.2`) and carry the rich invariant assertions; this harness is the
**shared cross-binding layer** — the same `fixtures.json`/`golden.json` a
non-Julia consumer reads — that asserts the documents reproduce a
binding-independent oracle.

## Layout

- `<model>/fixtures.json` — inputs (coords + projection parameters), shared across bindings
- `<model>/golden.json` — reference forward map + round-trip, an **independent closed-form oracle**
- `<model>/regenerate_golden.jl` — regenerates the golden from the closed form (no engine)
- `../../../test/test_reprojection_conformance.jl` — Julia (reference binding) consumer

## Models & invariants (RFC §8)

| Model              | Forward oracle                                   | Invariants asserted |
|--------------------|--------------------------------------------------|---------------------|
| `longlat`          | geographic WGS84 identity `x = lon − lon_0`      | forward == oracle (exact); inverse∘forward == identity; `lon_0` load-bearing |
| `lambert_conformal`| spherical LCC (Snyder 1987 = PROJ `+proj=lcc +R`)| forward == proj4 reference (WRF + NEI2016); origin → (0,0); inverse∘forward == identity to tol |

The golden's forward values are an **independent** implementation of the closed
form (`regenerate_golden.jl`), NOT the `.esm` AST: the LCC oracle is cross-checked
against the proj4-validated reference points from the co-located test before it is
frozen. Each binding then evaluates the `.esm` and must reproduce the oracle to the
fixture's declared tolerance.

## Cross-binding disposition (esd-47z.6)

Reprojection is a **point-wise scalar closed form** over the existing operator set
(`-`, `+`, `sin`, `cos`, `tan`, `atan`, `atan2`, `log`, `sqrt`, `sign`) — no
geometry kernel, no value-invention. Per RFC §8 the contract is
"`inverse ∘ forward ≈ identity` to tolerance … **cross-binding byte-identical
where the closed-form is exact**". The golden here is therefore the
binding-independent contract:

- **Reference binding (Julia):** evaluates the `.esm` through the ESS engine via
  the zero-IC surfacing idiom and asserts forward-vs-oracle + round-trip
  (`test/test_reprojection_conformance.jl`, this bead). **Green.**
- **Python / Rust / TypeScript:** the toolkit expression evaluators
  (`earthsci_toolkit` / `earthsci-toolkit`) already evaluate this operator set, so
  these consumers are mechanically portable against the SAME `fixtures.json` +
  `golden.json` (read the `.esm` forward/inverse observed ASTs, evaluate with the
  fixture bindings, compare to the oracle). They are tracked as **follow-up
  binding-port beads**, mirroring how the grid/DUO/structured cross-binding
  consumers were added incrementally after the Julia reference landed; no new
  operator or engine primitive is required.

Tolerances are declared per fixture. `longlat` is exact at the integer-valued
coords/`lon_0` (a 1e-12 floor guards cross-libm). `lambert_conformal` agrees with
the proj4 reference to `forward_relative = 1e-7` over the ~2000 km CONUS domain
(cross-libm rounding over a ~10-op chain) and round-trips to `1e-9`.

## Regenerating

```
julia tests/conformance/reprojection/longlat/regenerate_golden.jl
julia tests/conformance/reprojection/lambert_conformal/regenerate_golden.jl
```

Both activate a throwaway environment (only `JSON`); the goldens are pure closed
form, so regeneration needs neither the engine nor any binding.
