# earthsci_grids (Rust binding for EarthSciDiscretizations)

Rust implementation of the EarthSciDiscretizations grid accessors, conforming
to the cross-binding contract in [`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Build

This crate has a `path` dependency on the unpublished `earthsci-toolkit` crate
from [EarthSciSerialization](https://github.com/EarthSciML/EarthSciSerialization).
Check out ESS as a sibling of this repo:

```
parent/
├── EarthSciDiscretizations/   # this repo
└── EarthSciSerialization/     # https://github.com/EarthSciML/EarthSciSerialization
```

Then:

```bash
cargo build
cargo test
```

The `geo` and `proj` cargo features are declared placeholders for the
planned ecosystem integrations (`docs/GRIDS_API.md` §5.3); they pull in
the dependencies but gate no code yet, so enabling them currently adds
no functionality:

```bash
cargo build --features geo,proj
```

## Layout

- `src/lib.rs` — crate root, re-exports the public `Grid` trait and the
  per-family modules.
- `src/grids/` — per-family grid accessor modules (one file per family).
- `tests/` — integration tests.

See `docs/GRIDS_API.md` (repo root) for the normative API contract this
crate implements.
