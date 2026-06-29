#!/usr/bin/env python3
"""Generate conservative_regrid_overlap_join_computed.esm from the host-A_ij
conformance fixture conservative_regrid_overlap_join.esm.

The fixture is the worked example where the overlap matrix A_ij and the target
denominator dst_areas are supplied as HOST const_arrays (its conformance goldens
depend on that interface — see conservative_regrid_overlap_join_conformance_test.jl).
This sibling LANDS the declarative A_ij INTO the component: A_ij is COMPUTED from the
full source/target polygon rings via the ranged spherical intersect_polygon clip +
the same Van Oosterom-Strackee spherical-excess polygon_area FAQ the fixture already
uses for its representative pair — so no host A_ij/dst_areas. The broad-phase
bin-skolem join + sliver filter assembly is unchanged.

Reuses the fixture's proven VOS area expression verbatim (rewriting `clip[a,b]` ->
`clip[i,j,a,b]` so it ranges over candidate pairs) rather than re-deriving the
spherical-excess algebra. Output validates + runs through build_evaluator; verified
to reproduce the fixture's A_j/F_tgt on the same 3x3 spherical case.

Usage:  python3 gen_computed_overlap_join.py
"""
import copy
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "conservative_regrid_overlap_join.esm")
OUT = os.path.join(HERE, "conservative_regrid_overlap_join_computed.esm")


def rewrite_clip_indices(node):
    """Deep-copy `node`, prefixing every index(clip, a, b) -> index(clip, i, j, a, b)
    so the single-pair VOS area body ranges over the (i,j) candidate pairs."""
    if isinstance(node, dict):
        if node.get("op") == "index" and node.get("args") and node["args"][0] == "clip":
            rest = [rewrite_clip_indices(a) for a in node["args"][1:]]
            return {"op": "index", "args": ["clip", "i", "j"] + rest}
        return {k: rewrite_clip_indices(v) for k, v in node.items()}
    if isinstance(node, list):
        return [rewrite_clip_indices(x) for x in node]
    return node


def arrayop(oidx, ranges, expr):
    return {"op": "arrayop", "output_idx": list(oidx),
            "ranges": {k: {"from": v} for k, v in ranges},
            "args": [], "expr": expr}


def ix(*args):
    return {"op": "index", "args": list(args)}


def replace_index_var(node, name, new_name):
    """Deep-copy `node`, retargeting every index(`name`, ...) -> index(`new_name`, ...)
    and any bare `name` operand in an aggregate args list -> `new_name`."""
    if isinstance(node, dict):
        if node.get("op") == "index" and node.get("args") and node["args"][0] == name:
            return {"op": "index", "args": [new_name] + [replace_index_var(a, name, new_name)
                                                          for a in node["args"][1:]]}
        return {k: ([new_name if x == name else replace_index_var(x, name, new_name) for x in v]
                    if k == "args" and isinstance(v, list) else replace_index_var(v, name, new_name))
                for k, v in node.items()}
    if isinstance(node, list):
        return [replace_index_var(x, name, new_name) for x in node]
    return node


def main():
    doc = json.load(open(SRC))
    m = doc["models"]["ConservativeRegridOverlapJoin"]
    fv, mv = m["variables"], {e.get("_comment", ""): e for e in m["equations"]}

    # --- ranged spherical clip: clip[i,j,w,c] = intersect_polygon(src_poly[i], tgt_poly[j])[w,c]
    clip_expr = arrayop(
        ["i", "j", "w", "c"],
        [("i", "src_cells"), ("j", "tgt_cells"), ("w", "clip_ring"), ("c", "coord")],
        ix({"op": "intersect_polygon", "id": "overlap_clip", "manifold": "spherical",
            "args": [ix("src_poly", "i"), ix("tgt_poly", "j")]}, "w", "c"))

    # --- ranged A_ij: reuse the fixture's representative VOS area aggregate, with the
    #     clip refs lifted to (i,j), wrapped in an arrayop over the candidate pairs.
    vos_agg = rewrite_clip_indices(copy.deepcopy(fv["A_ij_rep"]["expression"]))
    A_ij_expr = arrayop(["i", "j"], [("i", "src_cells"), ("j", "tgt_cells")], vos_agg)

    # --- A_j_w: build-once target denominator Σ_i A_ij (plain sum; non-overlapping
    #     pairs contribute 0, so this equals the join-gated row-sum — proven shape).
    A_jw_expr = arrayop(["j"], [("j", "tgt_cells")],
                        {"op": "aggregate", "semiring": "sum_product", "output_idx": [],
                         "ranges": {"i": {"from": "src_cells"}}, "args": [],
                         "expr": ix("A_ij", "i", "j")})

    # --- Assembly: reuse the fixture's A_j / F_tgt EQUATIONS verbatim (the i,j
    #     double-range sum_product with the buffer-gated join + sliver filter, output
    #     j) — they reference A_ij by NAME, which is now the COMPUTED observed, so A_j
    #     needs no change. F_tgt: swap its host dst_areas denominator for the computed
    #     A_j_w observed (the only edit). Keeping the fixture's exact aggregate shape
    #     keeps the join-key resolution (tgt_bin -> range j) intact.
    eq_by_tag = {}
    for e in m["equations"]:
        c = e.get("_comment", "")
        if c.startswith("(3)"):
            eq_by_tag["A_j"] = copy.deepcopy(e)
        elif c.startswith("(4)+(5)"):
            eq_by_tag["F_tgt"] = copy.deepcopy(e)
    F_tgt_eq = eq_by_tag["F_tgt"]
    F_tgt_eq["rhs"] = replace_index_var(F_tgt_eq["rhs"], "dst_areas", "A_j_w")
    F_tgt_eq["_comment"] = ("(4)+(5) APPLY + NORMALIZE — F_tgt[j] = Σ_i A_ij·F_src[i] / A_j_w[j], "
                            "the COMPUTED target denominator (no host dst_areas). Same bin-skolem "
                            "buffer-gated join + sub-atol sliver filter as the fixture; only the "
                            "normalize divisor is now the computed A_j_w observed.")
    A_j_eq = eq_by_tag["A_j"]
    A_j_eq["_comment"] = ("(3) A_j[j] = Σ_i A_ij[i,j] — group-by-j row-sum over the bin-skolem "
                          "CANDIDATE pairs (join.on [[src_bin,tgt_bin]]) + sliver filter, now "
                          "reading the COMPUTED A_ij observed. Reproduces A_j_w (PoU by construction).")

    def P(*shape):
        return {"type": "parameter", "shape": list(shape)}

    def O(expr, *shape):
        return {"type": "observed", "shape": list(shape), "expression": expr}

    model = {
        "index_sets": copy.deepcopy(m["index_sets"]),
        "variables": {
            # binning representative coords (kept — the proven value-invention front-door)
            "src_lon": fv["src_lon"], "src_lat": fv["src_lat"],
            "tgt_lon": fv["tgt_lon"], "tgt_lat": fv["tgt_lat"],
            "dx": fv["dx"], "dy": fv["dy"], "atol": fv["atol"], "F_src": fv["F_src"],
            # full polygon rings (replace the representative src_poly_rep/tgt_poly_rep)
            "src_poly": P("src_cells", "cell_verts", "coord"),
            "tgt_poly": P("tgt_cells", "cell_verts", "coord"),
            # value-invention broad phase (kept verbatim)
            "src_bin": fv["src_bin"], "tgt_bin": fv["tgt_bin"], "pair_exists": fv["pair_exists"],
            # COMPUTED geometry (the landed declarative A_ij)
            "clip": O(clip_expr, "src_cells", "tgt_cells", "clip_ring", "coord"),
            "A_ij": O(A_ij_expr, "src_cells", "tgt_cells"),
            "A_j_w": O(A_jw_expr, "tgt_cells"),
            "A_j": fv["A_j"], "F_tgt": fv["F_tgt"],
        },
        "equations": [
            mv["(1a) OVERLAP PAIRS — bin each SOURCE cell to an integer lat-lon bin key. floor(lon/dx), floor(lat/dy) are integer bin coordinates (the existing `floor` op); skolem mints the content-addressed bin key. Per-element map (output i, no contraction)."],
            mv["(1b) OVERLAP PAIRS — bin each TARGET cell identically (same dx/dy quantization)."],
            None,  # candidate_set (1c) — filled below by id (its comment is long)
            A_j_eq,
            F_tgt_eq,
        ],
    }
    # candidate_set equation (1c): keep verbatim from the fixture
    for e in m["equations"]:
        if isinstance(e.get("rhs"), dict) and e["rhs"].get("id") == "candidate_set":
            model["equations"][2] = copy.deepcopy(e)
            break
    assert model["equations"][2] is not None, "candidate_set equation not found"

    md = copy.deepcopy(doc["metadata"])
    md["name"] = "conservative_regrid_overlap_join_computed"
    md["description"] = (
        "COMPUTED-A_ij conservative regridder (sibling of conservative_regrid_overlap_join.esm). "
        "The overlap matrix A_ij is COMPUTED declaratively from the full source/target polygon "
        "rings (src_poly[src_cells,cell_verts,coord], tgt_poly[tgt_cells,cell_verts,coord]) via the "
        "ranged spherical intersect_polygon clip + the Van Oosterom-Strackee spherical-excess "
        "polygon_area FAQ — NOT supplied as a host const_array. The target denominator A_j_w is the "
        "build-once row-sum Σ_i A_ij (no host dst_areas). The broad-phase bin-skolem floor->skolem->"
        "distinct->equijoin candidate set and the sub-atol sliver filter are unchanged from the "
        "fixture. This is the self-contained regridder the camp_fire.esm coupling needs: given only "
        "the cell geometries it produces the conservative F_tgt with no host-supplied overlaps. "
        "Generated by gen_computed_overlap_join.py from the conformance fixture; reproduces its "
        "A_j/F_tgt on the same spherical case.")
    md["tags"] = sorted(set(md.get("tags", []) + ["computed-aij", "self-contained"]))

    out = {"esm": doc["esm"], "metadata": md,
           "models": {"ConservativeRegridOverlapJoinComputed": model}}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print("wrote", OUT)
    print("  vars:", len(model["variables"]), " eqs:", len(model["equations"]),
          " index_sets:", list(model["index_sets"].keys()))


if __name__ == "__main__":
    main()
