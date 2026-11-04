#!/usr/bin/env python3
"""generate_catalog.py — build the generated docs content from the library files.

Walks grids/, regridding/, and reprojection/ and emits, under docs/content/:

  grids/<grid>/index.md      one page per grid: description, metaparameters,
                             index sets, geometry templates, stencils, and every
                             rule under grids/<grid>/rules/ rendered in place
                             (match pattern + discretization pretty-printed,
                             per-region boundary treatment, references, fixture
                             links, observed-order tables from committed
                             convergence goldens).
  regridding/_index.md       one section page listing entries with their
  reprojection/_index.md     templates pretty-printed.

and, under docs/static/plots/convergence/, a log-log convergence plot per
convergence case whose golden/errors.json exists (matplotlib; skipped with a
warning when matplotlib is absent).

Everything emitted here is generated ONLY in the CI docs job and never
committed (.gitignore pins docs/content/grids/, docs/content/regridding/,
docs/content/reprojection/, docs/static/plots/ — see AGENTS.md "Generated
files").

Expression display policy (AGENTS.md "single pathway"):
  * This script reads .esm JSON structurally and NEVER evaluates library math.
  * Every expression shown on a generated page is printed by the public
    `earthsci_ast.to_ascii` and nothing else, so the pages show EarthSciAST's
    ascii surface syntax — the text `earthsci_ast.parse_expression` reads —
    and a reader can paste any displayed expression back in. There is no
    local renderer and no fallback: earthsci_ast is a hard requirement (the CI
    docs job pip-installs it), and an expression `to_ascii` cannot print fails
    the build. Full precision always lives in the .esm sources the pages link.
  * Template names and parameter lists are shown as labels next to the
    printed body; region tables and stencil expansion are gone — a rule's
    `makearray` body is printed whole, imported stencils appearing as the
    template applications `name<param=value, ...>` they are.
  * Observed convergence orders are read from committed golden errors.json
    files (the same log-ratio arithmetic as scripts/check_convergence_order.py);
    no numbers are ever recomputed through any non-official pathway.

Usage:
  python3 docs/generate_catalog.py            generate into docs/
  python3 docs/generate_catalog.py --check    dry-run into a temp dir; exit
                                              nonzero if generation would fail
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path

DOCS = Path(__file__).resolve().parent
REPO = DOCS.parent
REPO_URL = "https://github.com/EarthSciML/earthscidiscretizations"

GRIDS = REPO / "grids"
REGRIDDING = REPO / "regridding"
REPROJECTION = REPO / "reprojection"
CONFORMANCE = REPO / "tests" / "conformance"

MINUS = "−"  # −
SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

warnings: list[str] = []
errors: list[str] = []


def warn(msg: str) -> None:
    warnings.append(msg)
    print(f"[warn] {msg}")


def fail(msg: str) -> None:
    errors.append(msg)
    print(f"[error] {msg}")



def load_json(path: Path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"{path.relative_to(REPO)}: unreadable JSON: {exc}")
        return None


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def blob(path: Path) -> str:
    return f"{REPO_URL}/blob/main/{rel(path)}"


def tree(path: Path) -> str:
    return f"{REPO_URL}/tree/main/{rel(path)}"


def tag_value(tags, prefix: str):
    for t in tags or []:
        if t.startswith(prefix):
            return t[len(prefix):]
    return None


def first_sentence(text: str) -> str:
    if not text:
        return ""
    for stop in (". ", ".\n"):
        i = text.find(stop)
        if i != -1:
            return text[: i + 1]
    return text


def yaml_str(value) -> str:
    """JSON is a subset of YAML — safe quoting for frontmatter scalars/lists."""
    return json.dumps(value, ensure_ascii=False)


def md_cell(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def md_prose(text: str) -> str:
    """Author prose from a `.esm` `metadata.description`, made safe to emit as Markdown.

    These descriptions are mathematics, not markup, and they are full of index
    subscripts. A subscript immediately followed by a parenthesised group —
    `kf[3/2](u[2]-u[1])/dx^2` — is *exactly* CommonMark inline-link syntax, so
    the renderer emits `<a href="u[2]-u[1]">3/2</a>` and the docs job's link
    checker then fails on a file that was never meant to be a link. Thirteen
    such phrases across three grid pages did precisely that.

    Escaping the closing bracket stops the link from forming and renders as the
    literal `]` the author wrote. Text inside a backtick code span is left
    alone: no link forms there, so escaping would only show a stray backslash.
    Nothing here can break a real link — citations live in
    `metadata.references`, and no `.esm` description in the library contains
    Markdown link syntax.
    """
    return _MD_CODE_OR_LINKISH.sub(
        lambda m: m.group("code") if m.group("code") else "\\](", text)


# Scanned left to right, alternation first: a backtick code span is consumed
# whole (so a `](` inside one is passed over), and only a `](` reached outside
# any code span is escaped. Splitting on code spans instead would mis-pair the
# backticks these descriptions already use around field names (`u`, `k`, `gi`)
# and escape the wrong halves.
_MD_CODE_OR_LINKISH = re.compile(r"(?P<code>`+[^`]*`+)|(?P<linkish>\]\()")


# ---------------------------------------------------------------------------
# Expression display: earthsci_ast.to_ascii, and nothing else.
# ---------------------------------------------------------------------------


def _load_to_ascii():
    try:
        # PUBLIC surface only, from the package namespace. In the Python
        # binding an `Expr` IS the wire dict, so raw .esm JSON goes in directly.
        from earthsci_ast import to_ascii  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "generate_catalog.py: earthsci_ast is required — every expression on the "
            "generated pages is printed by earthsci_ast.to_ascii. Install it with\n"
            '  pip install "git+https://github.com/EarthSciML/EarthSciAST.git'
            '#subdirectory=pkg/earthsci-ast-py"\n'
            f"({exc})"
        ) from exc
    return to_ascii


class AsciiPrinter:
    """Prints displayed expressions with the public `earthsci_ast.to_ascii`."""

    def __init__(self) -> None:
        self._to_ascii = _load_to_ascii()
        self.count = 0

    def __call__(self, node, where: str) -> str:
        try:
            text = self._to_ascii(node)
        except Exception as exc:  # reported, and fails the build
            fail(f"{where}: earthsci_ast.to_ascii failed: {exc}")
            return "?"
        self.count += 1
        return text


def code_span(text: str) -> str:
    """Inline Markdown code span outside a table (newlines folded)."""
    return "`" + text.replace("\n", " ") + "`"


# ---------------------------------------------------------------------------
# Library model: load grids, stencils, rules, cross-grid entries, fixtures.
# ---------------------------------------------------------------------------


def esm_files(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.iterdir() if p.suffix == ".esm" and p.is_file())


def sole_template(doc: dict, path: Path):
    """(name, template) for the single expression_templates entry, or None."""
    templates = doc.get("expression_templates") or {}
    if not templates:
        warn(f"{rel(path)}: no expression_templates")
        return None
    name = doc.get("metadata", {}).get("name") or next(iter(templates))
    if name not in templates:
        name = next(iter(templates))
    return name, templates[name]


# ---------------------------------------------------------------------------
# Consumer-supplied free names (the keyed-factor / free-name geometry contract).
#
# A grid's real-valued geometry is NOT baked into its .esm files: names like
# `x0`/`dx` (cartesian), `lon0_deg`/`dlon_deg` (latlon), or `areaCell`/`dvEdge`
# (mpas) appear as bare references in the geometry-template, stencil, and rule
# bodies and resolve in the CONSUMING model's scope at evaluation (esm-spec §9.6
# scalar-field substitution / the keyed-factor contract). They are neither
# metaparameters (load-time integers) nor template parameters (bound per apply
# site) nor loop variables — so we surface them structurally by walking the
# bodies and subtracting everything that IS bound. Presentation only; nothing
# here evaluates library math.
# ---------------------------------------------------------------------------

_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

# Dict keys whose string values are NOT expression name references (op names,
# axis/label metadata, template refs, reduction tags).
_NON_REF_KEYS = {
    "op", "wrt", "dim", "name", "semiring", "reduce", "manifold",
    "output_idx", "kind", "type", "units", "description", "priority",
}


def _walk_name_refs(node, bound: set, out: set) -> None:
    """Collect bare identifier references in expression positions not in `bound`."""
    if isinstance(node, str):
        if _IDENT.match(node) and node not in bound:
            out.add(node)
        return
    if isinstance(node, list):
        for item in node:
            _walk_name_refs(item, bound, out)
        return
    if not isinstance(node, dict):
        return
    local = bound
    if node.get("op") == "faq":
        # output_idx + range keys bind loop variables in this subtree.
        loop = set(node.get("output_idx") or []) | set((node.get("ranges") or {}).keys())
        local = bound | loop
    for key, value in node.items():
        if key in _NON_REF_KEYS:
            continue
        _walk_name_refs(value, local, out)


def free_names_for_grid(grid_doc: dict, sub_docs: list) -> list:
    """Consumer-supplied free names across a grid's geometry, stencils, and rules.

    `sub_docs` is the list of (path, doc) for the grid's stencils + rules.
    """
    bound: set = set()
    bound |= set((grid_doc.get("metaparameters") or {}).keys())
    index_sets = grid_doc.get("index_sets") or {}
    bound |= set(index_sets.keys())

    docs = [grid_doc] + [d for _, d in sub_docs]
    for doc in docs:
        templates = doc.get("expression_templates") or {}
        bound |= set(templates.keys())
        for tmpl in templates.values():
            bound |= set(tmpl.get("params") or [])

    # Ragged index sets name their keyed-factor arrays (offsets/values) — always
    # part of the consumer contract, even when they double as a template param
    # elsewhere (e.g. mpas nEdgesOnCell) or a body never spells them.
    keyed_factors: set = set()
    for spec in index_sets.values():
        if isinstance(spec, dict) and spec.get("kind") == "ragged":
            for k in ("offsets", "values"):
                v = spec.get(k)
                if isinstance(v, str) and _IDENT.match(v):
                    keyed_factors.add(v)

    body_refs: set = set()
    for doc in docs:
        for tmpl in (doc.get("expression_templates") or {}).values():
            body = tmpl.get("body")
            if body is not None:
                _walk_name_refs(body, bound, out=body_refs)
    body_free = {n for n in body_refs if n not in bound}
    return sorted(keyed_factors | body_free)


def index_conformance():
    """[{category, case, dir, manifest}] for every manifest under tests/conformance/."""
    cases = []
    if not CONFORMANCE.is_dir():
        return cases
    for cat_dir in sorted(p for p in CONFORMANCE.iterdir() if p.is_dir()):
        for case_dir in sorted(p for p in cat_dir.iterdir() if p.is_dir()):
            manifest_path = case_dir / "manifest.json"
            if not manifest_path.exists():
                continue
            manifest = load_json(manifest_path)
            if manifest is None:
                continue
            cases.append(
                {
                    "category": cat_dir.name,
                    "case": manifest.get("case", case_dir.name),
                    "dir": case_dir,
                    "manifest": manifest,
                }
            )
    return cases


def resolve_manifest_path(case_dir: Path, ref: str) -> Path | None:
    try:
        return (case_dir / ref).resolve()
    except OSError:
        return None


def problem_imports(problem_path: Path, cache: dict) -> set[Path]:
    """Resolved template-import refs (transitively shallow) of a problem file."""
    if problem_path in cache:
        return cache[problem_path]
    refs: set[Path] = set()
    doc = load_json(problem_path) if problem_path.exists() else None
    if doc:
        for model in (doc.get("models") or {}).values():
            for imp in model.get("expression_template_imports") or []:
                ref = imp.get("ref")
                if ref and not ref.startswith(("http:", "https:")):
                    try:
                        refs.add((problem_path.parent / ref).resolve())
                    except OSError:
                        pass
    cache[problem_path] = refs
    return refs


def fixtures_for_rule(rule_path: Path, cases, problem_cache) -> list[dict]:
    """Conformance cases that exercise this rule (by rule ref, ast dir name, or
    the problem file the case drives)."""
    resolved = rule_path.resolve()
    hits = []
    for case in cases:
        manifest, case_dir = case["manifest"], case["dir"]
        matched = False
        rule_ref = manifest.get("rule")
        if rule_ref and resolve_manifest_path(case_dir, rule_ref) == resolved:
            matched = True
        elif case["category"] == "ast" and case_dir.name == rule_path.stem:
            matched = True
        else:
            problem_ref = manifest.get("problem")
            if problem_ref:
                problem_path = resolve_manifest_path(case_dir, problem_ref)
                if problem_path and resolved in problem_imports(problem_path, problem_cache):
                    matched = True
        if matched:
            hits.append(case)
    return hits


def read_golden(case: dict):
    """Golden errors.json rows for a convergence case, or None (pending)."""
    manifest = case["manifest"]
    golden_path = case["dir"] / manifest.get("golden", "golden/errors.json")
    if not golden_path.exists():
        return None, golden_path
    golden = load_json(golden_path)
    if golden is None or "errors" not in golden:
        fail(f"{rel(golden_path)}: malformed convergence golden")
        return None, golden_path
    return sorted(golden["errors"], key=lambda r: r.get("n", 0)), golden_path


def observed_orders(ns, errs) -> list:
    """Same pairwise log-ratio as scripts/check_convergence_order.py."""
    out = [None]
    for k in range(1, len(errs)):
        if errs[k - 1] > 0 and errs[k] > 0 and ns[k] != ns[k - 1]:
            out.append(math.log(errs[k - 1] / errs[k]) / math.log(ns[k] / ns[k - 1]))
        else:
            out.append(None)
    return out


# ---------------------------------------------------------------------------
# Plots (matplotlib; degrade to a warning when absent).
# ---------------------------------------------------------------------------

_matplotlib_missing_warned = False


def write_convergence_plot(case: dict, rows, out_png: Path) -> bool:
    global _matplotlib_missing_warned
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        if not _matplotlib_missing_warned:
            warn("matplotlib not importable — convergence plots skipped (CI installs it)")
            _matplotlib_missing_warned = True
        return False

    manifest = case["manifest"]
    norms = manifest.get("norms") or [k for k in rows[0] if k != "n"]
    ns = [r["n"] for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    first_errs = None
    for norm in norms:
        errs = [r.get(norm) for r in rows]
        if any(e is None or e <= 0 for e in errs):
            continue
        if first_errs is None:
            first_errs = errs
        ax.loglog(ns, errs, marker="o", label=norm)
    p = manifest.get("expected_order")
    if p and first_errs:
        ref = [first_errs[-1] * (ns[-1] / n) ** p for n in ns]
        ax.loglog(ns, ref, "k--", linewidth=1, label=f"slope {MINUS}{p:g} reference")
    ax.set_xlabel("resolution n")
    ax.set_ylabel("error norm")
    ax.set_title(case["case"])
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)
    ax.legend(fontsize=8)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Page emission.
# ---------------------------------------------------------------------------


def math_block(text: str) -> str:
    return f'<div class="math-block">{html.escape(text)}</div>'


def chips(pairs) -> str:
    inner = " ".join(f"<code>{html.escape(str(p))}</code>" for p in pairs if p)
    return f'<p class="rule-meta">{inner}</p>'


def references_md(metadata: dict) -> list[str]:
    out = []
    refs = metadata.get("references") or []
    if not refs:
        return out
    out.append("**References**")
    out.append("")
    for ref in refs:
        citation = ref.get("citation", "")
        doi = ref.get("doi")
        line = f"- {md_cell(citation)}"
        if doi:
            # Percent-encode: DOIs legally contain <>() which break markdown links.
            from urllib.parse import quote

            line += f" [`doi:{doi}`](https://doi.org/{quote(doi, safe='/')})"
        out.append(line)
    out.append("")
    return out


def convergence_section(case: dict, plots_root: Path) -> list[str]:
    manifest = case["manifest"]
    name = case["case"]
    out = [f"#### Convergence — case [`{name}`]({tree(case['dir'])})", ""]
    rows, golden_path = read_golden(case)
    expected = manifest.get("expected_order")
    tolerance = manifest.get("order_tolerance")
    if rows is None:
        status = manifest.get("status", "pending-golden")
        out.append(
            f'<div class="callout callout-pending">Convergence golden not yet generated '
            f"(<code>status: {html.escape(str(status))}</code>). The observed-order table and "
            f"log-log plot will appear here automatically once "
            f"<code>{html.escape(rel(case['dir']))}/golden/errors.json</code> lands "
            f"(regenerated Julia-reference-only via <code>scripts/regenerate-goldens.sh</code>).</div>"
        )
        out.append("")
        if expected is not None:
            out.append(f"Expected order: **{expected:g}** (± {tolerance:g})." if tolerance is not None else f"Expected order: **{expected:g}**.")
            out.append("")
        return out

    norms = manifest.get("norms") or [k for k in rows[0] if k != "n"]
    ns = [r["n"] for r in rows]
    header = "| n |" + "".join(f" {n} | observed order |" for n in norms)
    sep = "|---|" + "---|---|" * len(norms)
    out += [
        f"Error norms read from the committed golden "
        f"[`{rel(golden_path)}`]({blob(golden_path)}) "
        f"(binding: `{manifest.get('reference_binding', 'julia')}`; nothing recomputed for display).",
        "",
        header,
        sep,
    ]
    orders = {norm: observed_orders(ns, [r.get(norm, 0) for r in rows]) for norm in norms}
    for k, row in enumerate(rows):
        cells = [str(row["n"])]
        for norm in norms:
            err = row.get(norm)
            cells.append(f"{err:.3e}" if isinstance(err, (int, float)) else "—")
            p = orders[norm][k]
            cells.append(f"{p:.2f}" if p is not None else "—")
        out.append("| " + " | ".join(cells) + " |")
    out.append("")
    if expected is not None:
        tol = f" (± {tolerance:g})" if tolerance is not None else ""
        out.append(f"Expected order: **{expected:g}**{tol}.")
        out.append("")
    plot_path = plots_root / "convergence" / f"{name}.png"
    if write_convergence_plot(case, rows, plot_path):
        out.append(f"![Convergence of {name} (log-log)](/plots/convergence/{name}.png)")
        out.append("")
    return out


def rule_section(
    rule_path: Path,
    doc: dict,
    printer: AsciiPrinter,
    conformance_cases,
    problem_cache,
    plots_root: Path,
) -> list[str]:
    metadata = doc.get("metadata", {})
    tags = metadata.get("tags") or []
    entry = sole_template(doc, rule_path)
    if entry is None:
        return []
    name, template = entry

    axes_tag = tag_value(tags, "axes:")

    out = [f"### `{name}`", ""]
    out.append(
        chips(
            [
                f"op:{tag_value(tags, 'op:')}" if tag_value(tags, "op:") else None,
                f"order:{tag_value(tags, 'order:')}" if tag_value(tags, "order:") else None,
                f"bc:{tag_value(tags, 'bc:')}" if tag_value(tags, "bc:") else None,
                f"axes:{axes_tag}" if axes_tag else None,
                f"stencil width {tag_value(tags, 'stencil_width:')}" if tag_value(tags, "stencil_width:") else None,
                f"priority {template['priority']}" if "priority" in template else None,
            ]
        )
    )
    out.append("")
    out.append(f"Source: [`{rel(rule_path)}`]({blob(rule_path)})")
    out.append("")

    # Match pattern.
    match = template.get("match")
    if match is not None:
        pretty = printer(match, f"{rel(rule_path)} match")
        priority_note = f" at priority {template['priority']}" if "priority" in template else ""
        out.append(f"**Rewrites** {code_span(pretty)}{priority_note} — match pattern:")
        out.append("")
        out.append("```json")
        out.append(json.dumps(match, ensure_ascii=False))
        out.append("```")
        out.append("")
    else:
        warn(f"{rel(rule_path)}: rule template has no match pattern")

    # Match-scoping `where` shape constraint (esm-spec §9.6.1): the rule fires
    # only on a bare field of the named shape, which is also the mechanism that
    # lets the rule be imported multiple times (renamed) without collision.
    where = template.get("where")
    if where:
        clauses = []
        for var, spec in where.items():
            shape = spec.get("shape") if isinstance(spec, dict) else None
            if isinstance(shape, list):
                axes_s = ", ".join(str(a) for a in shape)
                clauses.append(f"<code>{html.escape(str(var))}</code> is a bare field shaped <code>[{html.escape(axes_s)}]</code>")
            else:
                clauses.append(f"<code>{html.escape(str(var))}</code>: <code>{html.escape(json.dumps(spec, ensure_ascii=False))}</code>")
        out.append(
            '<div class="callout callout-scope"><strong>Match scope</strong> '
            "(esm-spec §9.6.1 <code>where</code>): fires only when " + "; ".join(clauses) + ". "
            "Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape "
            "follow the renamed axis together, so the rule can be imported more than once "
            "(each instance scoped to its own grid) without a first-declared-wins collision. "
            "A consumer differentiating a compound inline expression must bind it to a "
            "declared shaped observed first.</div>"
        )
        out.append("")

    # Discretization: the whole body, printed by earthsci_ast.to_ascii.
    body = template.get("body")
    if body is not None:
        if isinstance(body, dict) and body.get("op") == "makearray":
            out.append(
                "**Discretization** — one `makearray` tiling the output axes; boundary "
                "conditions live in the face regions of this same rule (esm-spec §9.6.8; "
                "later regions overwrite earlier, §4.3.2). Imported stencils appear as "
                "template applications `name<param=value, ...>`, defined in the "
                "Stencils section of this page:"
            )
        else:
            out.append("**Definition:**")
        out.append("")
        out.append(math_block(printer(body, f"{rel(rule_path)} body")))
        out.append("")

    description = metadata.get("description")
    if description:
        out.append(md_prose(description))
        out.append("")

    out += references_md(metadata)

    # Conformance fixtures on GitHub.
    hits = fixtures_for_rule(rule_path, conformance_cases, problem_cache)
    out.append("**Conformance fixtures**")
    out.append("")
    if hits:
        for case in hits:
            status = case["manifest"].get("status")
            status_note = f" — status `{status}`" if status else ""
            out.append(
                f"- {case['category']}: [`{case['case']}`]({tree(case['dir'])}){status_note}"
            )
    else:
        out.append("- none yet")
    out.append("")

    for case in hits:
        if case["category"] == "convergence":
            out += convergence_section(case, plots_root)

    return out


def params_line(params) -> str:
    """Label line naming a template's parameters (shown beside its printed body)."""
    if not params:
        return "No parameters."
    return "Parameters: " + ", ".join(code_span(str(p)) for p in params)


def stencil_section(stencil_path: Path, doc: dict, printer: AsciiPrinter) -> list[str]:
    metadata = doc.get("metadata", {})
    tags = metadata.get("tags") or []
    entry = sole_template(doc, stencil_path)
    if entry is None:
        return []
    name, template = entry
    out = [f"### `{name}`", ""]
    out.append(
        chips(
            [
                f"op:{tag_value(tags, 'op:')}" if tag_value(tags, "op:") else None,
                f"order:{tag_value(tags, 'order:')}" if tag_value(tags, "order:") else None,
                f"stencil width {tag_value(tags, 'stencil_width:')}" if tag_value(tags, "stencil_width:") else None,
                "interior-only (match-less)",
            ]
        )
    )
    out.append("")
    out.append(f"Source: [`{rel(stencil_path)}`]({blob(stencil_path)})")
    out.append("")
    out.append(params_line(template.get("params")))
    out.append("")
    body = template.get("body")
    if body is not None:
        out.append(math_block(printer(body, f"{rel(stencil_path)} body")))
        out.append("")
    description = metadata.get("description")
    if description:
        out.append(md_prose(description))
        out.append("")
    out += references_md(metadata)
    return out


def grid_page(
    grid_dir: Path,
    printer: AsciiPrinter,
    conformance_cases,
    problem_cache,
    plots_root: Path,
) -> str | None:
    grid_path = grid_dir / "grid.esm"
    doc = load_json(grid_path)
    if doc is None:
        return None
    metadata = doc.get("metadata", {})
    grid_name = metadata.get("name", grid_dir.name)
    tags = metadata.get("tags") or []
    family = tag_value(tags, "family:")

    stencil_docs = [(p, load_json(p)) for p in esm_files(grid_dir / "stencils")]
    rule_docs = [(p, load_json(p)) for p in esm_files(grid_dir / "rules")]
    stencil_docs = [(p, d) for p, d in stencil_docs if d is not None]
    rule_docs = [(p, d) for p, d in rule_docs if d is not None]

    all_tags = list(tags)
    ops, bcs = set(), set()
    for _, d in stencil_docs + rule_docs:
        dtags = d.get("metadata", {}).get("tags") or []
        all_tags += [t for t in dtags if t not in all_tags]
        op = tag_value(dtags, "op:")
        bc = tag_value(dtags, "bc:")
        if op:
            ops.add(op)
        if bc:
            bcs.add(bc)
    families = sorted(
        {family} | {tag_value(d.get("metadata", {}).get("tags") or [], "family:") for _, d in stencil_docs + rule_docs} - {None}
    ) if family else sorted(
        {tag_value(d.get("metadata", {}).get("tags") or [], "family:") for _, d in stencil_docs + rule_docs} - {None}
    )

    lines = [
        "---",
        f"title: {yaml_str(grid_name)}",
        f"description: {yaml_str(first_sentence(metadata.get('description', '')))}",
        "esd_kind: \"grid\"",
        f"families: {yaml_str(families)}",
        f"ops: {yaml_str(sorted(ops))}",
        f"bcs: {yaml_str(sorted(bcs))}",
        f"tags: {yaml_str(all_tags)}",
        f"source: {yaml_str(rel(grid_path))}",
        "---",
        "",
    ]
    if metadata.get("description"):
        lines.append(md_prose(metadata["description"]))
        lines.append("")

    metaparameters = doc.get("metaparameters") or {}
    if metaparameters:
        lines += ["## Metaparameters", "", "| Name | Type | Default | Description |", "|---|---|---|---|"]
        for pname, spec in metaparameters.items():
            default = spec.get("default", "—")
            lines.append(
                f"| `{pname}` | {spec.get('type', '')} | `{default}` | {md_cell(md_prose(spec.get('description', '')))} |"
            )
        lines += [
            "",
            "All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the "
            "import edge (or through the loader API) and the same files serve every resolution.",
            "",
        ]

    free_names = free_names_for_grid(doc, stencil_docs + rule_docs)
    if free_names:
        recipe = {
            "cartesian": "The geometry is consumer-supplied so the same files serve any "
            "**domain extent**, not just any resolution. For a domain [a, b] the consuming "
            "model defines `x0 = a` (a parameter) and `dx = (b − a)/N` (an observed whose "
            "expression divides by the metaparameter name `N`, so a loader-API rebinding of "
            "`N` — a convergence sweep — keeps `dx` consistent). A model that instead closes "
            "`N` at the import edge spells the matching literal (`{op: /, args: [1, 8]}` for "
            "N = 8 on the unit interval). See the extent proof "
            "`problems/heat_1d_zero_grad_nonunit.esm` (x ∈ [−1.5, 2.5], observed order 2.00).",
            "latlon": "Region-generic: the same files serve a **global** or a **regional** "
            "grid depending on the free-name values the consumer supplies. GLOBAL recipe — "
            "`lon0_deg = 0`, `dlon_deg = 360/NLON`, `lat0_deg = −90`, `dlat_deg = "
            "180/(NLAT − 1)` (pole-to-pole; the zonal circle closes, so the periodic-lon rule "
            "applies). REGIONAL recipe — any other origin/spacing (e.g. `lon0_deg = −130`, "
            "`dlon_deg = 0.25` over CONUS); the circle does not close, so use the "
            "zero-gradient lon rule. Spell the spacings as observeds dividing by the "
            "metaparameter names (`NLON`/`NLAT`) so a convergence sweep stays consistent. A "
            "further free name `R_sphere` (sphere radius, m) enters only the consumer's "
            "physical-units metric composition (per-degree derivatives → per-metre; the "
            "`coslat` template supplies cos(lat)).",
            "mpas": "The mesh geometry arrives as **data**, not formula: these are the "
            "keyed-factor arrays a consuming model must expose under exactly these bare names, "
            "normally as observed aliases of a mounted `grids/mpas/mesh/` subsystem's `const` "
            "data (see `problems/divergence_mpas.esm`). The 2026-07 regridding wave widened "
            "the contract with the vertex-geometry factors so cell polygon rings derive "
            "in-document via `mpas_cell_rings_deg`.",
        }.get(
            family,
            "These names are supplied by the consuming model's scope at evaluation (the "
            "keyed-factor / free-name geometry contract); they are neither metaparameters "
            "(load-time integers) nor template parameters (bound per apply site).",
        )
        chips_html = " ".join(f"<code>{html.escape(n)}</code>" for n in free_names)
        lines += [
            "## Consumer-supplied free names",
            "",
            f'<p class="rule-meta">{chips_html}</p>',
            "",
            recipe,
            "",
        ]

    index_sets = doc.get("index_sets") or {}
    if index_sets:
        lines += ["## Index sets", "", "| Name | Kind | Size |", "|---|---|---|"]
        for iname, spec in index_sets.items():
            if spec.get("kind") == "ragged":
                of = spec.get("of") or []
                detail = f"per {', '.join(of)}" if of else "ragged"
                extras = ", ".join(
                    f"{k}: {spec[k]}" for k in ("offsets", "values") if k in spec
                )
                size = f"{detail} ({extras})" if extras else detail
            else:
                size = str(spec.get("size", ""))
            lines.append(f"| `{iname}` | {spec.get('kind', '')} | `{md_cell(size)}` |")
        lines.append("")

    templates = doc.get("expression_templates") or {}
    if templates:
        lines += [
            "## Geometry templates",
            "",
            "| Template | Parameters | Definition | Description |",
            "|---|---|---|---|",
        ]
        for tname, tmpl in templates.items():
            params = ", ".join(f"`{p}`" for p in tmpl.get("params") or []) or "—"
            definition = printer(tmpl.get("body"), f"{rel(grid_path)} template {tname}")
            lines.append(
                f"| `{tname}` | {params} | `{md_cell(definition)}` "
                f"| {md_cell(md_prose(tmpl.get('description', '')))} |"
            )
        lines.append("")

    if stencil_docs:
        lines += [
            "## Stencils",
            "",
            "Interior-only, match-less templates importing `grid.esm`; the rules below wrap "
            "them with boundary-condition face regions into complete rewrite rules.",
            "",
        ]
        for path, sdoc in stencil_docs:
            lines += stencil_section(path, sdoc, printer)

    lines += ["## Rules", ""]
    if rule_docs:
        lines += [
            "Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary "
            "conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching "
            "derivative in your model lowers through it.",
            "",
        ]
        for path, rdoc in rule_docs:
            lines += rule_section(
                path, rdoc, printer, conformance_cases, problem_cache, plots_root
            )
    else:
        lines += ["No rules published for this grid yet.", ""]

    return "\n".join(lines) + "\n"


def cross_grid_section_page(directory: Path, title: str, kind: str, printer: AsciiPrinter) -> str:
    """Section page (_index.md) for regridding/ or reprojection/."""
    tag_prefix = {"regrid": "method:", "reproject": "crs:"}[kind]
    tag_label = {"regrid": "method", "reproject": "crs"}[kind]
    intro = {
        "regrid": "Conservative and interpolating remap expressions (coupling transforms). "
        "Each entry is a pure ESM template library — declarative expressions, no kernels.",
        "reproject": "Coordinate-transform template fragments (longitude-latitude, Lambert "
        "conformal, …) for use inside regridding and coupling expressions.",
    }[kind]

    capability = {
        "regrid": [
            "The conservative overlap regridder is **end-to-end declarative**: a consumer "
            "constructs source and target cell rings, derives the broad-phase bin keys, and "
            "runs the candidate-gated overlap → normalize → apply, all from `.esm` templates "
            "with no host glue.",
            "",
            "- **Cell-ring constructors** turn a grid spec into the `[cells, verts, coord]` "
            "vertex rings the regridder consumes — `cartesian_cell_rings.esm` from a uniform "
            "cartesian spec, and `grids/mpas/grid.esm`'s `mpas_cell_rings_deg` from MPAS mesh "
            "vertex data (the widened contract).",
            "- **Broad phase in-library**: representative binning coordinates derived from the "
            "rings themselves plus integer `skolem` bin keys, so the candidate set is "
            "byte-identical across bindings.",
            "- **Gated == dense** value-identically under an admissible binning (the gated "
            "denominator broad-phases like the numerator, so partition-of-unity stays exact).",
            "- **One file, every manifold**: `planar` or `spherical` is a per-invocation "
            "template parameter. Spherical clip runs on Julia + Rust today; the Python "
            "spherical case is gated on the optional `spherely` dependency (planar cases run "
            "on Julia + Python + Rust), and the two rewrite-only ports (Go, TypeScript) are "
            "scope-excluded — recorded in each case manifest, never shimmed.",
            "",
            "End-to-end exemplars: `tests/conformance/regridding/"
            "cartesian_rings_regrid_gated_3x3_to_2x2` (grid spec → rings → broad phase → "
            "gated overlap → regrid) and `mpas_l0_to_octants_sphere` (MPAS mesh-derived rings "
            "on the sphere).",
        ],
        "reproject": [
            "These are forward/inverse **scalar** templates over the evaluable core — but they "
            "also apply **in-model over coordinate arrays**: invoke them inside a `faq` node "
            "whose bindings mix indexed reads, literals, and parameter references, and every "
            "`apply_expression_template` inlines (esm-spec §9.6.2 Option A) into a closed "
            "scalar AST at each cell. The `tests/conformance/reprojection/lcc_grid_roundtrip` "
            "case pins exactly this faq-mapped lowering: byte-identical expanded AST "
            "across all five bindings, and a working round-trip simulation on Julia + Python "
            "+ Rust.",
        ],
    }[kind]

    entries = esm_files(directory)
    lines = [
        "---",
        f"title: {yaml_str(title)}",
        f"description: {yaml_str(intro)}",
        "---",
        "",
        intro,
        "",
    ] + capability + [""]
    if not directory.is_dir() or not entries:
        lines += [
            f'<div class="callout callout-pending">No <code>{directory.name}/</code> entries are '
            "published yet — this section fills in automatically as library files land "
            f'(see <a href="{tree(REPO / directory.name) if directory.is_dir() else REPO_URL}">the repository</a>).</div>',
            "",
        ]
        return "\n".join(lines) + "\n"

    for path in entries:
        doc = load_json(path)
        if doc is None:
            continue
        metadata = doc.get("metadata", {})
        name = metadata.get("name", path.stem)
        tags = metadata.get("tags") or []
        tagv = tag_value(tags, tag_prefix)
        lines += [f"### `{name}`", ""]
        lines.append(chips([f"{tag_label}:{tagv}" if tagv else None]))
        lines.append("")
        lines.append(f"Source: [`{rel(path)}`]({blob(path)})")
        lines.append("")
        if metadata.get("description"):
            lines.append(md_prose(metadata["description"]))
            lines.append("")
        templates = doc.get("expression_templates") or {}
        if templates:
            lines.append(f"<details><summary>Templates ({len(templates)})</summary>")
            lines.append("")
            for tname, tmpl in templates.items():
                lines.append(f"**`{tname}`**")
                if tmpl.get("description"):
                    lines.append("")
                    lines.append(md_prose(tmpl["description"]))
                lines.append("")
                lines.append(params_line(tmpl.get("params")))
                lines.append("")
                body = tmpl.get("body")
                if body is not None:
                    lines.append(math_block(printer(body, f"{rel(path)} template {tname}")))
                    lines.append("")
            lines.append("</details>")
            lines.append("")
        lines += references_md(metadata)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------


def clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def generate(content_root: Path, static_root: Path) -> None:
    printer = AsciiPrinter()

    conformance_cases = index_conformance()
    problem_cache: dict = {}
    plots_root = static_root / "plots"

    grids_out = content_root / "grids"
    clean_dir(grids_out)
    clean_dir(plots_root)

    (grids_out / "_index.md").write_text(
        "---\n"
        'title: "Grids"\n'
        'description: "One page per grid: index sets, geometry, stencils, and the complete '
        'boundary-condition-carrying rewrite rules published for it."\n'
        "---\n\n"
        "Everything specific to a grid lives under `grids/<grid>/` — `grid.esm` (index sets "
        "+ geometry, metaparameterized), `stencils/` (interior operators), `rules/` (stencil "
        "+ boundary conditions = the complete rewrite rule). The same scheme takes a "
        "different form on different grids; each page below documents one grid and every "
        "rule available on it.\n",
        encoding="utf-8",
    )

    n_grids = 0
    if GRIDS.is_dir():
        for grid_dir in sorted(p for p in GRIDS.iterdir() if p.is_dir()):
            if not (grid_dir / "grid.esm").exists():
                warn(f"{rel(grid_dir)}: no grid.esm — skipped (possibly still landing)")
                continue
            page = grid_page(grid_dir, printer, conformance_cases, problem_cache, plots_root)
            if page is None:
                continue
            out_dir = grids_out / grid_dir.name
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "index.md").write_text(page, encoding="utf-8")
            n_grids += 1
    else:
        warn("grids/ directory missing — no grid pages generated")

    for directory, title, kind in (
        (REGRIDDING, "Regridding", "regrid"),
        (REPROJECTION, "Reprojection", "reproject"),
    ):
        out_dir = content_root / directory.name
        clean_dir(out_dir)
        (out_dir / "_index.md").write_text(
            cross_grid_section_page(directory, title, kind, printer), encoding="utf-8"
        )

    print(f"[info] generated {n_grids} grid page(s) into {grids_out}")
    print(f"[info] earthsci_ast.to_ascii printed {printer.count} expression(s)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="generate into a throwaway directory; exit nonzero if generation would fail",
    )
    args = parser.parse_args()

    if args.check:
        with tempfile.TemporaryDirectory(prefix="esd-docs-check-") as tmp:
            generate(Path(tmp) / "content", Path(tmp) / "static")
    else:
        generate(DOCS / "content", DOCS / "static")

    if warnings:
        print(f"[summary] {len(warnings)} warning(s)")
    if errors:
        print(f"[summary] {len(errors)} error(s) — generation FAILED")
        return 1
    print("[summary] OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
