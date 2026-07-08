# `flux_match` — DECLARATIVE GAP (esd-6g4.7 / G7)

**Verdict: the `flux_match` interface modifier CANNOT currently be expressed as a
declarative rewrite rule that fires through the ESS engine. Reported per the
bead's DECLARATIVE-OR-FAIL directive — NOT implemented imperatively.**

`zero_gradient` and `interface` (this directory's sibling rules) DO fire — see
`zero_gradient_bc.json` and `interface_bc.json` with their `fixtures/rewrite/`
goldens. `flux_match` is the one of the three named kinds that is blocked, and
the block is precise and external (an ESS matcher gap), not a modeling choice.

## What `flux_match` is

Per `esm-spec.md` §11.5 (Boundary Condition Types) and `esm-schema.json`,
`flux_match` is **not a boundary-condition `kind`**. The closed set of BC kinds
is `{constant, zero_gradient, periodic, dirichlet, neumann, robin, interface}`.
`flux_match` is a **boolean field that modifies an `interface` BC**:

> `flux_match` | boolean | For `kind='interface'`: if true, enforce normal-flux
> continuity (∂u/∂n = ∂v/∂n) in addition to value continuity at the interface.
> Default false.

So a `flux_match` BC is really `kind: "interface"` plus `flux_match: true`.

## Why it cannot fire through the current engine

The ESS BC kind/side matcher (G8, landed in `ess-bps`) works by having
`_discretize_bc!` (EarthSciSerialization `src/discretize.jl`) build a **synthetic
`bc` wrapper** from the BC spec and run the rule engine on it. The wrapper
encodes ONLY these fields:

| BC spec field      | Wrapper encoding                          | Matched by            |
|--------------------|-------------------------------------------|-----------------------|
| `variable`         | first `args` entry                        | structural arg match  |
| `kind`             | `fn` sibling field                        | `_match_sibling_name` |
| `side`             | `dim` sibling field                       | `_match_sibling_name` |
| `coupled_variable` | second `args` entry                       | structural arg match  |
| `value`            | trailing `args` entry (when present)      | structural arg match  |

`flux_match` is **never read** by `_discretize_bc!` and **never placed on the
`bc` wrapper**. (Confirmed on ESS `origin/main`: the only Julia reference to
`flux_match` is the schema; `grep flux_match packages/EarthSciSerialization.jl/src`
returns nothing. The toolkit bindings — Python/Go/TS/Rust — parse and serialize
the field, but the Julia discretizer does not surface it to pattern matching.)

Therefore there is **no field a rule pattern could match on** to distinguish an
`interface` BC with `flux_match: true` from one with `flux_match: false`. Both
lower to the identical wrapper `bc[fn=interface, dim=side](var, coupled)`, so any
rule that matched the flux-matching variant would fire identically on the plain
variant — the modifier's semantics (the additional ∂u/∂n = ∂v/∂n constraint)
cannot be realized. A rule keyed on a hypothetical `fn: "flux_match"` is also
dead: `flux_match` is not a valid BC `kind`, so the wrapper's `fn` is never that
string for any schema-valid model.

Independently, the math reinforces the gap: meaningful flux matching across an
interface with *differing* transport coefficients needs those coefficients
(a harmonic-mean / coefficient-weighted ghost), and those values are likewise
absent from the matchable wrapper.

## What would unblock it (ESS-side, not ESD-side)

A future ESS extension to `_discretize_bc!` could propagate the modifier into the
matchable wrapper, e.g. one of:

1. Emit a **distinct synthetic `fn`** when `flux_match` is true (e.g.
   `fn = "interface_flux_match"`), so a dedicated declarative rule can match it;
   or
2. Carry `flux_match` (and, for the coefficient-weighted case, the relevant
   transport coefficients) as **additional matchable `bc` fields/args**.

Either is a generic ESS engine change (mirrors how `coupled_variable` was added
as the second arg for `interface` in `ess-bps`). Once it lands, the
`flux_match` rule + a `fixtures/rewrite/` golden can be authored here exactly
like `interface_bc.json` — no ESD-side imperative operator. File against ESS.

## Scope note

Per DECLARATIVE-OR-FAIL, no imperative/hand-coded `flux_match` operator was
written. The two feasible kinds (`zero_gradient`, `interface`) ship as firing
declarative rules + fixtures; this verdict is the deliverable for the third.
