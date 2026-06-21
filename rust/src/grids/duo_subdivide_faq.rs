//! DUO icosahedral mesh construction with every vertex COORDINATE routed through
//! the ESS front-door evaluator (`eval_coeff`), preserving the EXACT combinatorial
//! ORDER of the imperative `subdivide_icosahedron` (`grids::duo`).
//!
//! Rust mirror of `src/subdivide_faq.jl` (D3, esd-ohd / W1). Declarative
//! companions: `discretizations/grids/duo/faq/subdivide_fold.esm` and
//! `subdivide_refine.esm`. The seed `vert_raw / sqrt(seed_sq)` normalize and the
//! midpoint `(a+b) / sqrt(Σ (a+b)²)` average are the geometry FAQs this file
//! evaluates through `eval_coeff` (`crate::rule_eval` — the single arithmetic
//! pathway; no shadow evaluator, per AGENTS.md / GRIDS_API §4.3).
//!
//! ## Front-door coordinates, imperative encounter-order numbering
//!
//! Only the **coordinates** go through the front door: every vertex value (the 12
//! seed vertices and every refined midpoint) is computed by `eval_coeff` on the
//! SAME ASTs the Julia binding drives. The **numbering**, however, is the
//! imperative encounter order — 12 base vertices in table order, then midpoints
//! minted in `HashMap`-cache ENCOUNTER order as faces are iterated — NOT the
//! canonical sorted/skolem/rank order serialized in `subdivide_refine.esm`. The
//! production path must keep the imperative order because downstream
//! Voronoi-dual conformance asserts the per-vertex dual arrays equal the
//! imperative dual EXACTLY, and the dual is indexed by primal vertex.
//!
//! Because the Rust `eval_coeff` path and the imperative `duo.rs` path call the
//! same Rust `f64` libm (`sqrt`, the left-fold `+`/`*`), the result is
//! BYTE-IDENTICAL and same-order as `subdivide_icosahedron(level)`,
//! coordinate-for-coordinate.

use std::collections::HashMap;

use serde_json::{json, Value};

use crate::{eval_coeff, Bindings};

/// Raw golden-ratio icosahedron vertices, pre-normalization. Transcribed to match
/// `icosahedron_vertices_f64`'s `raw` literal table (`grids::duo`) row-for-row;
/// phi = (1 + sqrt(5))/2.
fn faq_vert_raw() -> [[f64; 3]; 12] {
    let phi = (1.0f64 + 5.0f64.sqrt()) / 2.0;
    [
        [0.0, 1.0, phi],
        [0.0, -1.0, phi],
        [0.0, 1.0, -phi],
        [0.0, -1.0, -phi],
        [1.0, phi, 0.0],
        [-1.0, phi, 0.0],
        [1.0, -phi, 0.0],
        [-1.0, -phi, 0.0],
        [phi, 0.0, 1.0],
        [phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
    ]
}

// --- Seed normalize-to-sphere FAQ AST ----------------------------------------
// coord_d = r$d / sqrt(r1*r1 + r2*r2 + r3*r3). Squares are products (`*`), NOT
// `^`, and the 3-arg `+` left-folds — the exact float-op sequence of the
// imperative `n = sqrt(v0*v0 + v1*v1 + v2*v2); v[d]/n`.
fn rsq(name: &str) -> Value {
    json!({"op": "*", "args": [name, name]})
}

fn seed_norm() -> Value {
    json!({"op": "sqrt", "args": [{"op": "+", "args": [rsq("r1"), rsq("r2"), rsq("r3")]}]})
}

fn seed_comp() -> [Value; 3] {
    [
        json!({"op": "/", "args": ["r1", seed_norm()]}),
        json!({"op": "/", "args": ["r2", seed_norm()]}),
        json!({"op": "/", "args": ["r3", seed_norm()]}),
    ]
}

// --- Midpoint normalize-to-sphere FAQ AST ------------------------------------
// coord_d = (a$d + b$d) / sqrt(Σ_d (a$d + b$d)*(a$d + b$d)). Squares are products
// so the float ops match the imperative `midpoint` bit-for-bit.
fn sum_d(d: usize) -> Value {
    json!({"op": "+", "args": [format!("a{d}"), format!("b{d}")]})
}

fn sq_d(d: usize) -> Value {
    json!({"op": "*", "args": [sum_d(d), sum_d(d)]})
}

fn mid_norm() -> Value {
    json!({"op": "sqrt", "args": [{"op": "+", "args": [sq_d(1), sq_d(2), sq_d(3)]}]})
}

fn mid_comp() -> [Value; 3] {
    [
        json!({"op": "/", "args": [sum_d(1), mid_norm()]}),
        json!({"op": "/", "args": [sum_d(2), mid_norm()]}),
        json!({"op": "/", "args": [sum_d(3), mid_norm()]}),
    ]
}

/// The 12 base icosahedron vertices on the unit sphere, in `faq_vert_raw` table
/// order (= `icosahedron_vertices_f64` order), each computed via the seed
/// normalize FAQ (`eval_coeff`). Bit-identical to `icosahedron_vertices_f64()`
/// for `f64`.
fn faq_seed_vertices(seed: &[Value; 3]) -> Vec<[f64; 3]> {
    let mut verts: Vec<[f64; 3]> = Vec::with_capacity(12);
    for r in faq_vert_raw().iter() {
        let mut b = Bindings::new();
        b.insert("r1".to_string(), r[0]);
        b.insert("r2".to_string(), r[1]);
        b.insert("r3".to_string(), r[2]);
        verts.push([
            eval_coeff(&seed[0], &b).expect("duo seed faq: eval_coeff"),
            eval_coeff(&seed[1], &b).expect("duo seed faq: eval_coeff"),
            eval_coeff(&seed[2], &b).expect("duo seed faq: eval_coeff"),
        ]);
    }
    verts
}

/// Encounter-order edge-midpoint minting, mirroring the imperative `midpoint`
/// EXACTLY in cache keying (sorted pair `a < b`) and numbering (next index on
/// first encounter), but with the midpoint COORDINATE computed via the midpoint
/// normalize FAQ (`eval_coeff`). Bit-identical to `midpoint` for `f64`.
fn faq_midpoint(
    verts: &mut Vec<[f64; 3]>,
    cache: &mut HashMap<(u32, u32), u32>,
    a: u32,
    b: u32,
    mid: &[Value; 3],
) -> u32 {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&i) = cache.get(&key) {
        return i;
    }
    let va = verts[a as usize];
    let vb = verts[b as usize];
    let mut bnd = Bindings::new();
    bnd.insert("a1".to_string(), va[0]);
    bnd.insert("a2".to_string(), va[1]);
    bnd.insert("a3".to_string(), va[2]);
    bnd.insert("b1".to_string(), vb[0]);
    bnd.insert("b2".to_string(), vb[1]);
    bnd.insert("b3".to_string(), vb[2]);
    verts.push([
        eval_coeff(&mid[0], &bnd).expect("duo midpoint faq: eval_coeff"),
        eval_coeff(&mid[1], &bnd).expect("duo midpoint faq: eval_coeff"),
        eval_coeff(&mid[2], &bnd).expect("duo midpoint faq: eval_coeff"),
    ]);
    let idx = (verts.len() - 1) as u32;
    cache.insert(key, idx);
    idx
}

/// Build the level-`level` DUO icosahedral mesh on the UNIT sphere, returning
/// `(verts, faces)` — a drop-in replacement for the imperative
/// `subdivide_icosahedron(level)` (`grids::duo`). Every coordinate is computed
/// through the ESS front-door evaluator `eval_coeff`, but the vertex numbering and
/// face order are the imperative ENCOUNTER order — NOT the canonical skolem/rank
/// order. Byte-identical (`f64`) and same-order as `subdivide_icosahedron(level)`.
pub(crate) fn duo_subdivide_faq(level: u32) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let seed = seed_comp();
    let mid = mid_comp();

    let mut verts = faq_seed_vertices(&seed);
    // Pure-integer structural face table reused from `grids::duo` (no floats, so
    // the reuse is exact — no transcription risk).
    let mut faces: Vec<[u32; 3]> = super::duo::icosahedron_faces_0based().to_vec();

    // Refine: replicate `subdivide_icosahedron`'s loop EXACTLY — fresh cache per
    // level, iterate faces in order, mint ab/bc/ca via the encounter-order
    // midpoint FAQ, push the 4 children in the imperative order.
    for _ in 0..level {
        let mut cache: HashMap<(u32, u32), u32> = HashMap::new();
        let mut new_faces: Vec<[u32; 3]> = Vec::with_capacity(faces.len() * 4);
        for f in &faces {
            let a = f[0];
            let b = f[1];
            let c = f[2];
            let ab = faq_midpoint(&mut verts, &mut cache, a, b, &mid);
            let bc = faq_midpoint(&mut verts, &mut cache, b, c, &mid);
            let ca = faq_midpoint(&mut verts, &mut cache, c, a, &mid);
            new_faces.push([a, ab, ca]);
            new_faces.push([b, bc, ab]);
            new_faces.push([c, ca, bc]);
            new_faces.push([ab, bc, ca]);
        }
        faces = new_faces;
    }
    (verts, faces)
}
