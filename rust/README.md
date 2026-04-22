# earthsci_grids (Rust binding for EarthSciDiscretizations)

Rust implementation of the EarthSciDiscretizations grid accessors, conforming
to the cross-binding contract in [`docs/GRIDS_API.md`](../docs/GRIDS_API.md).

## Build

```bash
cargo build
cargo test
```

Optional ecosystem features (per `docs/GRIDS_API.md` §5.3):

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
