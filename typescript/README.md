# @earthsci/grids (TypeScript binding for EarthSciDiscretizations)

TypeScript implementation of the EarthSciDiscretizations grid accessors,
conforming to the cross-binding contract in
[`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Install

This package depends on the unpublished `earthsci-toolkit` package from
[EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization)
via a `file:` reference. Check out ESS as a sibling to this repo, then install:

```
parent/
├── EarthSciDiscretizations/   # this repo
└── EarthSciSerialization/     # https://github.com/EarthSciML/EarthSciSerialization
```

```bash
# from parent/
git clone https://github.com/EarthSciML/EarthSciSerialization.git
(cd EarthSciSerialization/packages/earthsci-toolkit && npm install && npm run build)

# then in this directory:
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
