# Reduction cross-binding conformance harness

Verifies the declarative dimension-reduction models (`reduction/<model>.esm`) are
green **through the ESS engine** and pins the binding-neutral fixtures + golden
that every binding consumes.

Scope: bead `esd-khv` (the campfire-e2e C4 design spike, 2026-06-26). The C4
pipeline lands loader fields on the Camp Fire target grid via reproject + regrid;
the 3-D ERA5 fields are first collapsed to the ground surface with `lev=min`
**before** the horizontal regrid. This family owns that vertical reduction. The
per-model **Julia evaluator** test ships co-located with the `.esm`
(`reduction/lev_min_surface_reduce_conformance_test.jl`) and carries the rich
end-to-end assertions; this harness is the **shared cross-binding layer**.

## Layout

- `<model>/fixtures.json` — inputs (3-D field + vertical coordinate), shared across bindings
- `<model>/golden.json` — reference surface field + the argmin level, an **independent analytic oracle**
- `<model>/regenerate_golden.jl` — regenerates the golden from the closed form (no engine)
- `../../../test/test_reduction_conformance.jl` — Julia (reference binding) consumer

## Models & invariants

| Model     | Reduction                                                        | Cross-binding mode |
|-----------|-----------------------------------------------------------------|--------------------|
| `lev_min` | `F_surf[x,y] = F_3d[x,y, argmin_lev lev_coord]` (value:min slice)| **binding-independent** (exact argmin gather) |

## Cross-binding disposition (esd-khv)

### `lev_min` — binding-independent, exact

The ground_surface interface's `dimension_mapping.constraints.lev = {value:min}`
expressed PURELY over existing aggregate vocab — an outer `sum_product` that
carries the spatial output `(x,y)` and contracts `lev`, whose body keeps a level
iff its coordinate equals the **inline** `min_sum` reduction `MIN_k lev_coord[k]`:

```
F_surf[x,y] = SUM_lev ifelse(lev_coord[lev] == MIN_k lev_coord[k], F_3d[x,y,lev], 0)
```

NO new engine op (the C4 declarative-or-fail gate, "lev=min reduces to
`aggregate`", **passes**); NO value-invention front-door (unlike `point`); NO
geometry clip (unlike `conservative`). This is the same tier as `bspline`
(regridding): a **closed-form oracle every binding can reproduce**. It is also
**exact (0-ULP)** rather than tolerance-based: `min` returns one of the
`lev_coord` elements bit-for-bit, so the engine's `==` indicator selects the
argmin level byte-exactly and the surfaced value is a verbatim `F_3d` gather. The
`absolute` tolerance in `fixtures.json` is a guard band, not a consumed budget.

What IS load-bearing in the shared fixture: the minimum is placed at an
**interior** level (lev 2 of 3), so a binding that shortcuts `lev=min` to a
first- or last-slice produces a different field and fails the golden — the argmin
search is genuinely exercised.

**Precondition:** a UNIQUE minimum (a strictly-monotone vertical coordinate —
physically the case for pressure / height / model-level). A tie would surface the
sum of the tied levels; interfaces declare a single `value:min` level.

The Julia reference consumer is **green**; Python/Rust/TS consumers are
mechanically portable **follow-ups** (the oracle is a plain argmin gather, and an
engine-backed binding drives the same `.esm`).

## Regenerating

```
julia tests/conformance/reduction/lev_min/regenerate_golden.jl   # analytic (JSON only)
```

The golden regenerates from a throwaway `JSON`-only environment (no engine).
