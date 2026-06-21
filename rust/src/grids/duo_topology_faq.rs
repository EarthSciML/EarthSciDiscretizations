//! DUO primal triangular-mesh topology via the landed ESS value-invention engine
//! (`earthsci_toolkit` relational primitives).
//!
//! Rust mirror of `src/topology_faq.jl` (D1a). Declarative companion:
//! `discretizations/grids/duo/rules/primal_topology.esm` expresses this chain as a
//! value-invention FAQ document (RFC `semiring-faq-unified-ir` §7.3: enumerate
//! edges from the face→vertex relation, `skolem` on the sorted vertex-pair key,
//! `distinct` dedup, `rank` dense ids, then the inversion join for
//! `cell_neighbors` and the group-by for `vertex_faces`).
//!
//! This module is the thin **evaluation bridge**: it shapes a triangular `faces`
//! table into the relational inputs and routes every dedup / rank / join through
//! ESS's `relational` primitives. ESD hosts no shadow relational engine — the
//! determinism contract (and the cross-binding byte-identity of the output) lives
//! entirely in `earthsci_toolkit` (CONFORMANCE_SPEC §5.5; AGENTS.md
//! single-pathway invariant). The only ESD-side logic is *authoring* which vertex
//! pairs are a triangle's edges.
//!
//! Edge numbering is the ESS canonical (sorted) order produced by `distinct` —
//! the same order the imperative `build_connectivity` already emits (it
//! `sort_unstable()`s the edge keys), so the edge array is byte-identical too.
//! `cell_neighbors` (cell×slot, `NO_NEIGHBOR` sentinel) and `vertex_faces` are
//! byte-identical regardless of edge numbering.

use earthsci_toolkit::{distinct, equijoin, rank, skolem_edge, Key};

use super::duo::NO_NEIGHBOR;
use crate::{GridError, Result};

/// Materialized DUO primal unstructured connectivity (matches the `DuoGrid`
/// field conventions).
pub(crate) struct PrimalTopology {
    /// (Ne) sorted vertex pairs `(a < b)`, in canonical (sorted) order.
    pub edges: Vec<[u32; 2]>,
    /// (Nc) cell sharing edge `k` of face `c`; `NO_NEIGHBOR` for a boundary edge
    /// (none on a closed icosahedral mesh).
    pub cell_neighbors: Vec<[u32; 3]>,
    /// Per-vertex sorted incident-face ids.
    pub vertex_faces: Vec<Vec<u32>>,
}

/// The undirected canonical edge key for a vertex pair (sorted `(min, max)`).
fn edge_key(a: u32, b: u32) -> Key {
    skolem_edge(Key::Int(a as i64), Key::Int(b as i64))
}

/// The three undirected edges of a triangle, indexed by the local slot `k`
/// OPPOSITE vertex `k` (k=0→(v1,v2), k=1→(v2,v0), k=2→(v0,v1)) — the slot
/// convention shared by `cell_neighbors` and matched to `build_connectivity`.
fn face_edges(f: &[u32; 3]) -> [(u32, u32); 3] {
    [(f[1], f[2]), (f[2], f[0]), (f[0], f[1])]
}

fn key_int(k: &Key) -> i64 {
    match k {
        Key::Int(i) => *i,
        _ => unreachable!("relational topology: expected Key::Int component"),
    }
}

fn tuple_int(k: &Key, i: usize) -> i64 {
    match k {
        Key::Tuple(v) => key_int(&v[i]),
        _ => unreachable!("relational topology: expected Key::Tuple row"),
    }
}

/// equijoin key projection: the dense edge id (component 0 of the incidence row).
fn join_on_edge_id(k: &Key) -> Key {
    match k {
        Key::Tuple(v) => v[0].clone(),
        _ => unreachable!("relational topology: expected Key::Tuple incidence row"),
    }
}

/// Materialize DUO primal connectivity from a `(Nc)` triangular `faces` table
/// (0-based vertex ids) by evaluating the value-invention FAQ through the landed
/// ESS relational engine. `distinct` (dedup), `rank` (dense ids), and `equijoin`
/// (connectivity inversion) are computed by `earthsci_toolkit`, so the result is
/// byte-identical to the Julia / Python bindings' engine by ESS's cross-binding
/// determinism contract.
pub(crate) fn primal_topology_faq(faces: &[[u32; 3]], nv: usize) -> Result<PrimalTopology> {
    let nc = faces.len();

    // --- (1) edges: skolem on the sorted vertex pair + distinct + rank ----------
    let mut edge_keys: Vec<Key> = Vec::with_capacity(3 * nc);
    for f in faces {
        for (a, b) in face_edges(f) {
            edge_keys.push(edge_key(a, b));
        }
    }
    let edge_order = distinct(&edge_keys); // sorted unique (min, max) pairs
    let ranking = rank(&edge_keys); // 0-based dense ids (Rust native)
    let mut edges: Vec<[u32; 2]> = Vec::with_capacity(edge_order.len());
    for k in &edge_order {
        edges.push([tuple_int(k, 0) as u32, tuple_int(k, 1) as u32]);
    }

    // --- (2) cell_neighbors: the OTHER cell across each edge --------------------
    // Inversion join (ESS equijoin): self-join the (dense edge id, cell, slot)
    // incidence on the edge id; each pair (c1, c2) with c1 != c2 sharing an edge
    // is the neighbour across that slot.
    let mut incidence: Vec<Key> = Vec::with_capacity(3 * nc);
    for (c, f) in faces.iter().enumerate() {
        for (k, (a, b)) in face_edges(f).into_iter().enumerate() {
            let eid = ranking.id[&edge_key(a, b)];
            incidence.push(Key::Tuple(vec![
                Key::Int(eid),
                Key::Int(c as i64),
                Key::Int(k as i64),
            ]));
        }
    }
    let joined = equijoin(&incidence, &incidence, join_on_edge_id, join_on_edge_id);
    let mut neighbors = vec![NO_NEIGHBOR; 3 * nc];
    for (l, r) in &joined {
        let c1 = tuple_int(l, 1) as usize;
        let k1 = tuple_int(l, 2) as usize;
        let c2 = tuple_int(r, 1) as u32;
        if c1 as u32 == c2 {
            continue;
        }
        let slot = 3 * c1 + k1;
        if neighbors[slot] != NO_NEIGHBOR && neighbors[slot] != c2 {
            return Err(GridError::SchemaViolation(format!(
                "duo: non-manifold edge — cell {c1} slot {k1} shares an edge with >2 cells"
            )));
        }
        neighbors[slot] = c2;
    }
    let cell_neighbors: Vec<[u32; 3]> = (0..nc)
        .map(|c| [neighbors[3 * c], neighbors[3 * c + 1], neighbors[3 * c + 2]])
        .collect();

    // --- (3) vertex_faces: faces incident on each vertex, sorted ----------------
    // Group-by over the face→vertex relation: `distinct` over (vertex, face) rows
    // sorts by (vertex, face), so each vertex's faces land contiguous and ascending.
    let mut vf_rows: Vec<Key> = Vec::with_capacity(3 * nc);
    for (c, f) in faces.iter().enumerate() {
        for &v in f.iter() {
            vf_rows.push(Key::Tuple(vec![Key::Int(v as i64), Key::Int(c as i64)]));
        }
    }
    let vf_sorted = distinct(&vf_rows);
    let mut vertex_faces: Vec<Vec<u32>> = vec![Vec::new(); nv];
    for row in &vf_sorted {
        let v = tuple_int(row, 0) as usize;
        let f = tuple_int(row, 1) as u32;
        vertex_faces[v].push(f);
    }

    Ok(PrimalTopology {
        edges,
        cell_neighbors,
        vertex_faces,
    })
}
