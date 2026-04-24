//! Cross-language conformance test for the vertical grid family.
//!
//! Every binding (python, typescript, julia, rust) must emit byte-identical
//! canonical `.esm` content for the same vertical-grid options — per
//! `docs/GRIDS_API.md` §3.5 and §4.1. The vertical family has no
//! transcendental math, so the §4.2 ULP fallback is not expected to fire;
//! `tolerance.relative` in `fixtures.json` is 0.0 (strict byte equality).
//!
//! The committed `discretizations/grids/vertical/*.esm` files are both the
//! ESD fixture corpus and the conformance golden — duplicating them under a
//! local `golden/` would leave two sources of truth for the same
//! declarative config. Mirrors `python/tests/test_vertical_conformance.py`,
//! `typescript/tests/vertical.conformance.test.ts`, and
//! `test/test_vertical_conformance.jl`.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::vertical::{self, VerticalCoordinate};
use earthsci_grids::Grid;
use serde_json::Value;

fn harness_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR is .../rust; harness lives at ../tests/conformance/...
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/vertical")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
}

fn coordinate_from(name: &str) -> VerticalCoordinate {
    VerticalCoordinate::from_name(name).unwrap_or_else(|| panic!("unknown coordinate {name}"))
}

fn f64_array(v: &Value, field: &str) -> Vec<f64> {
    v.as_array()
        .unwrap_or_else(|| panic!("{field} must be an array"))
        .iter()
        .map(|x| {
            x.as_f64()
                .unwrap_or_else(|| panic!("{field} element must be numeric"))
        })
        .collect()
}

/// Build the vertical grid for a fixture's declarative `opts`. Keeps the
/// same accepted key set as the python / julia / typescript sibling
/// runners.
fn build_grid(opts: &Value) -> earthsci_grids::vertical::VerticalGrid {
    let coordinate = coordinate_from(
        opts["coordinate"]
            .as_str()
            .expect("fixture opts.coordinate must be a string"),
    );
    let mut b = vertical::builder().coordinate(coordinate);
    if let Some(nz) = opts.get("nz") {
        b = b.nz(nz.as_u64().expect("nz must be integer") as usize);
    }
    if let Some(levels) = opts.get("levels") {
        b = b.levels(f64_array(levels, "levels"));
    }
    if let Some(ak) = opts.get("ak") {
        b = b.ak(f64_array(ak, "ak"));
    }
    if let Some(bk) = opts.get("bk") {
        b = b.bk(f64_array(bk, "bk"));
    }
    if let Some(p0) = opts.get("p0") {
        b = b.p0(p0.as_f64().expect("p0 must be numeric"));
    }
    if let Some(t) = opts.get("transition") {
        b = b.transition(t.as_f64().expect("transition must be numeric"));
    }
    if let Some(ghosts) = opts.get("ghosts") {
        b = b.ghosts(ghosts.as_u64().expect("ghosts must be integer") as usize);
    }
    b.build().expect("vertical grid build failed")
}

/// Render a finite `f64` with Python `repr`-compatible formatting. Rust's
/// `f64::to_string` emits `"1"` for `1.0`; Python emits `"1.0"`. Append the
/// trailing `.0` for integer-valued floats so the canonical bytes match
/// across bindings. No scientific-notation floats appear in the vertical
/// fixture set, so no additional normalisation is required (the julia and
/// typescript canonicalizers make the same assumption).
fn fmt_f64(x: f64) -> String {
    assert!(
        x.is_finite(),
        "vertical conformance canonicalization: non-finite number {x}"
    );
    let s = x.to_string();
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{s}.0")
    }
}

fn fmt_number(n: &serde_json::Number) -> String {
    if n.is_f64() {
        fmt_f64(n.as_f64().expect("f64 number"))
    } else {
        // Integers round-trip via Number's Display.
        n.to_string()
    }
}

/// Serialize a `Value` with object keys sorted recursively and 2-space
/// indent, matching Python's `json.dumps(..., sort_keys=True, indent=2)`.
fn serialize(v: &Value, indent: usize) -> String {
    match v {
        Value::Null => "null".to_string(),
        Value::Bool(b) => if *b { "true" } else { "false" }.to_string(),
        Value::Number(n) => fmt_number(n),
        Value::String(s) => serde_json::to_string(s).expect("json string"),
        Value::Array(xs) => {
            if xs.is_empty() {
                return "[]".to_string();
            }
            let inner = " ".repeat(indent + 2);
            let parts: Vec<String> = xs
                .iter()
                .map(|x| format!("{inner}{}", serialize(x, indent + 2)))
                .collect();
            format!("[\n{}\n{}]", parts.join(",\n"), " ".repeat(indent))
        }
        Value::Object(map) => {
            if map.is_empty() {
                return "{}".to_string();
            }
            let sorted: BTreeMap<&String, &Value> = map.iter().collect();
            let inner = " ".repeat(indent + 2);
            let parts: Vec<String> = sorted
                .iter()
                .map(|(k, val)| {
                    format!(
                        "{inner}{}: {}",
                        serde_json::to_string(k).expect("json string"),
                        serialize(val, indent + 2)
                    )
                })
                .collect();
            format!("{{\n{}\n{}}}", parts.join(",\n"), " ".repeat(indent))
        }
    }
}

/// Canonical Python-compatible JSON for the emitted Rust `Value`. Produces
/// exactly what `json.dumps(doc, sort_keys=True, indent=2) + "\n"` would
/// emit in Python for the same document. Values come from `to_esm()` on
/// a freshly-built grid, so every `f64` is the exact IEEE-754 result —
/// there is no round-trip through `serde_json`'s float parser, which on
/// long-digit inputs like `0.9683012089999999` can round to the neighbour
/// f64.
fn canonicalize_emitted(doc: &Value) -> String {
    serialize(doc, 0) + "\n"
}

/// Remove lines whose trimmed content starts with `"<key>":` for any `key`
/// in `strip_fields`. Cheap, robust, and — crucially — avoids round-
/// tripping the golden through `serde_json`'s f64 parser, whose default
/// (non-arbitrary-precision) path can land 1 ULP off exact values like the
/// eta-fixture's synthesized sigma. Per `fixtures.json`'s
/// `provenance_strip_fields` the stripped keys appear on their own lines
/// inside `provenance` at a deeper indent than any other sibling-key
/// occurrence of the same name; the check is therefore safe in the current
/// corpus. If a future fixture puts one of these names in a non-stripped
/// position, the line-based filter would over-strip and the equality check
/// would fail loudly, surfacing the conflict rather than silently passing.
fn strip_provenance_lines(text: &str, strip_fields: &[String]) -> String {
    let prefixes: Vec<String> = strip_fields.iter().map(|k| format!("\"{k}\":")).collect();
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        let trimmed = line.trim_start();
        if prefixes.iter().any(|p| trimmed.starts_with(p)) {
            continue;
        }
        out.push_str(line);
    }
    out
}

#[test]
fn vertical_matches_committed_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"]
        .as_f64()
        .expect("tolerance.relative must be numeric");
    assert_eq!(
        rel_tol, 0.0,
        "Vertical family is declarative-only; ULP fallback not expected. \
         If {rel_tol} is intentional, update this assertion and the other \
         runners together (docs/GRIDS_API.md §4.2 requires review across all \
         bindings)."
    );

    let strip_fields: Vec<String> = spec["provenance_strip_fields"]
        .as_array()
        .expect("provenance_strip_fields must be an array")
        .iter()
        .map(|v| v.as_str().expect("strip field must be string").to_string())
        .collect();
    let golden_dir = hdir.join(
        spec["golden_dir"]
            .as_str()
            .expect("golden_dir must be string"),
    );

    for fixture in spec["fixtures"].as_array().expect("fixtures array") {
        let name = fixture["name"].as_str().expect("fixture name");
        let grid = build_grid(&fixture["opts"]);
        let emitted_full = canonicalize_emitted(&grid.to_esm());
        let emitted = strip_provenance_lines(&emitted_full, &strip_fields);

        let golden_path = golden_dir.join(format!("{name}.esm"));
        let golden_text = fs::read_to_string(&golden_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", golden_path.display()));
        let expected = strip_provenance_lines(&golden_text, &strip_fields);

        assert_eq!(
            emitted, expected,
            "rust vertical binding disagrees with committed golden \
             ({name}.esm) after stripping provenance.{:?}",
            strip_fields
        );
    }
}

#[test]
fn every_committed_fixture_is_in_corpus() {
    // If a new fixture is committed under discretizations/grids/vertical/
    // without being added to fixtures.json, a binding could silently diverge
    // on it. This check is the tripwire, mirroring the python and julia
    // runners.
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let golden_dir = hdir.join(
        spec["golden_dir"]
            .as_str()
            .expect("golden_dir must be string"),
    );

    let mut committed: Vec<String> = fs::read_dir(&golden_dir)
        .expect("golden_dir readable")
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let p = e.path();
            if p.extension().and_then(|s| s.to_str()) == Some("esm") {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string())
            } else {
                None
            }
        })
        .collect();
    committed.sort();

    let mut corpus: Vec<String> = spec["fixtures"]
        .as_array()
        .expect("fixtures array")
        .iter()
        .map(|f| f["name"].as_str().expect("fixture name").to_string())
        .collect();
    corpus.sort();

    assert_eq!(
        committed, corpus,
        "Committed vertical fixtures diverged from the conformance corpus. \
         committed={committed:?}, corpus={corpus:?}. Add the new fixture to \
         tests/conformance/grids/vertical/fixtures.json or remove the orphan \
         .esm file."
    );
}
