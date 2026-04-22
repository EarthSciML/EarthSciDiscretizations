# @earthsci/grids (TypeScript binding for EarthSciDiscretizations)

TypeScript implementation of the EarthSciDiscretizations grid accessors,
conforming to the cross-binding contract in
[`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Install

```bash
npm install
```

## Layout

- `src/index.ts` — package entry; exposes `grids` namespace and `earthsci.grids` alias.
- `src/grids/` — per-family grid accessor modules (one file per family).
- `tests/` — vitest suite.

## Scripts

```bash
npm run typecheck    # tsc --noEmit
npm run build        # compile to dist/
npm test             # vitest run
npm run lint         # eslint
```

See `docs/GRIDS_API.md` (repo root) for the normative API contract this
package implements.
