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

Pretty-printing policy (AGENTS.md "single pathway"):
  * This script reads .esm JSON structurally and NEVER evaluates library math.
  * Float constants display at <= 12 significant digits ("%.12g", trailing
    zeros stripped) in BOTH display paths: that is the official
    earthsci_toolkit.display convention (_format_number), and the structural
    fallback here normalizes to the same so generated pages are byte-identical
    regardless of which renderer handled an expression. Presentation-layer
    only — full precision always lives in the .esm sources the pages link.
  * Math display prefers the OFFICIAL ESS display path
    (earthsci_toolkit.display.to_unicode), used in-process when the toolkit is
    importable (the CI docs job pip-installs it) or through
    docs/_render_bridge.py in the EarthSciSerialization toolkit venv when a
    local checkout is available.
  * Node types the official renderer does not display usefully — index
    arithmetic, aggregate, makearray, apply_expression_template, nested D — are
    rendered structurally here (compact math-ish text such as
    "(f[i+1] − 2·f[i] + f[i−1]) / dx²"). That is presentation, not evaluation.
  * With no official display available at all, everything falls back to the
    structural renderer and a warning is printed.
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
import os
import shutil
import subprocess
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
DOT = "·"  # ·
PARTIAL = "∂"  # ∂
ELEM = "∈"  # ∈
SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

warnings: list[str] = []
errors: list[str] = []


def warn(msg: str) -> None:
    warnings.append(msg)
    print(f"[warn] {msg}")


def fail(msg: str) -> None:
    errors.append(msg)
    print(f"[error] {msg}")


def sup(n) -> str:
    return str(n).translate(SUPERSCRIPT)


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


# ---------------------------------------------------------------------------
# Official ESS display path (in-process or via the toolkit venv bridge).
# ---------------------------------------------------------------------------


class OfficialDisplay:
    """Best-effort handle on earthsci_toolkit's to_unicode()."""

    def __init__(self) -> None:
        self.mode: str | None = None
        self._proc = None
        self._inproc = None
        try:
            from earthsci_toolkit.display import to_unicode  # type: ignore
            from earthsci_toolkit.parse import _parse_expression  # type: ignore

            self._inproc = (to_unicode, _parse_expression)
            self.mode = "in-process"
            return
        except Exception:
            pass
        self._try_bridge()

    def _try_bridge(self) -> None:
        ess_root = Path(os.environ.get("ESS_ROOT") or (REPO.parent / "EarthSciSerialization"))
        toolkit = ess_root / "packages" / "earthsci_toolkit"
        venv_py = toolkit / ".venv" / "bin" / "python"
        if not venv_py.exists():
            return
        env = dict(os.environ)
        env["PYTHONPATH"] = str(toolkit / "src") + os.pathsep + env.get("PYTHONPATH", "")
        try:
            proc = subprocess.Popen(
                [str(venv_py), str(DOCS / "_render_bridge.py")],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                env=env,
            )
            handshake = json.loads(proc.stdout.readline() or "{}")
            if handshake.get("ok"):
                self._proc = proc
                self.mode = f"ess-venv bridge ({venv_py})"
            else:
                proc.kill()
        except Exception:
            self._proc = None

    @property
    def available(self) -> bool:
        return self.mode is not None

    def render(self, expr) -> str | None:
        if self._inproc is not None:
            to_unicode, parse = self._inproc
            try:
                return to_unicode(parse(expr))
            except Exception:
                return None
        if self._proc is not None:
            try:
                self._proc.stdin.write(json.dumps(expr) + "\n")
                self._proc.stdin.flush()
                out = json.loads(self._proc.stdout.readline() or "{}")
                return out["text"] if out.get("ok") else None
            except Exception:
                return None
        return None

    def close(self) -> None:
        if self._proc is not None:
            try:
                self._proc.stdin.close()
                self._proc.wait(timeout=5)
            except Exception:
                self._proc.kill()


# ---------------------------------------------------------------------------
# Structural unicode renderer (presentation only — no evaluation).
# ---------------------------------------------------------------------------

# Ops the official display path cannot render usefully for library files
# (index(f, 2) instead of f[2]; apply_expression_template() losing the name;
# aggregate(f) hiding the kernel). These always render structurally.
STRUCTURAL_ONLY_OPS = {"index", "aggregate", "makearray", "apply_expression_template"}


def official_safe(node) -> bool:
    """True when the subtree contains no node the official renderer mangles."""
    if isinstance(node, (int, float, str)) or node is None:
        return True
    if isinstance(node, dict):
        op = node.get("op")
        if op in STRUCTURAL_ONLY_OPS:
            return False
        if op == "D":
            arg = (node.get("args") or [None])[0]
            if isinstance(arg, dict) and arg.get("op") == "D":
                return False  # official prints ∂∂f/∂x/∂x for nested D
        for key in ("args",):
            for a in node.get(key) or []:
                if not official_safe(a):
                    return False
        return True
    return False


class MathRenderer:
    SUM, TERM, UNARY, POW, ATOM = 40, 50, 60, 70, 100

    def __init__(self, official: OfficialDisplay) -> None:
        self.official = official
        self.official_uses = 0
        self.structural_fallback_ops: set[str] = set()

    # -- public entry points ------------------------------------------------

    def render_top(self, node) -> str:
        """Render a whole expression, preferring the official display path."""
        if self.official.available and official_safe(node):
            text = self.official.render(node)
            if text is not None:
                self.official_uses += 1
                return text
        return self.render(node, 0)

    REDUCTION_SYMBOL = {"sum_product": "Σ", "prod": "Π", "min_plus": "min", "max_plus": "max"}

    def render_aggregate(self, node) -> tuple[str, str, list[str]]:
        """-> (kernel expression, range description, output indices).

        Indices in `ranges` but not in `output_idx` are contracted (reduced by
        the semiring); they render under the reduction symbol, together with
        any `filter` condition.
        """
        output_idx = list(node.get("output_idx") or [])
        ranges = node.get("ranges") or {}
        expr = self.render_top(node.get("expr"))
        cond = self.render_top(node["filter"]) if node.get("filter") else None
        out_ranges = {k: v for k, v in ranges.items() if k in output_idx}
        contracted = {k: v for k, v in ranges.items() if k not in output_idx}
        if contracted:
            sym = self.REDUCTION_SYMBOL.get(node.get("semiring", "sum_product"), node.get("semiring", "Σ"))
            under = self.render_ranges(contracted)
            if cond:
                under += f" | {cond}"
            expr = f"{sym}_{{{under}}} {expr}"
            cond = None
        ranges_str = self.render_ranges(out_ranges)
        if cond:
            ranges_str += f" where {cond}"
        return expr, ranges_str, output_idx

    def render_ranges(self, ranges: dict) -> str:
        parts = []
        for name, spec in ranges.items():
            if isinstance(spec, dict) and "from" in spec:
                of = spec.get("of")
                suffix = f"({', '.join(of)})" if isinstance(of, list) and of else ""
                parts.append(f"{name} {ELEM} {spec['from']}{suffix}")
            elif isinstance(spec, (list, tuple)) and len(spec) == 2:
                lo, hi = (self.render(b, 0, compact=True) for b in spec)
                parts.append(f"{name} {ELEM} [{lo}, {hi}]")
            else:
                parts.append(f"{name} {ELEM} ?")
        return ", ".join(parts)

    def render_region(self, region, axes: list[str]) -> tuple[str, str]:
        """One makearray region -> (label, role) with role interior|face."""
        parts, face = [], False
        for k, bounds in enumerate(region):
            axis = axes[k] if k < len(axes) else f"d{k + 1}"
            lo, hi = (self.render(b, 0, compact=True) for b in bounds)
            if lo == hi:
                parts.append(f"{axis} = {lo}")
                face = True
            else:
                parts.append(f"{axis} {ELEM} [{lo}, {hi}]")
        return ", ".join(parts), ("boundary face" if face else "interior")

    # -- recursive structural rendering --------------------------------------

    def render(self, node, prec: int, compact: bool = False) -> str:
        if isinstance(node, bool):
            return str(node).lower()
        if isinstance(node, int):
            return self._number(node, prec)
        if isinstance(node, float):
            return self._number(node, prec)
        if isinstance(node, str):
            return node
        if not isinstance(node, dict):
            return str(node)

        op = node.get("op")
        args = node.get("args") or []

        if op == "+":
            return self._sum(args, prec, compact)
        if op == "-":
            if len(args) == 1:
                s = MINUS + self.render(args[0], self.UNARY, compact)
                return self._wrap(s, self.SUM, prec)
            left = self.render(args[0], self.SUM, compact)
            right = self.render(args[1], self.SUM + 2, compact)
            sep = MINUS if compact else f" {MINUS} "
            return self._wrap(f"{left}{sep}{right}", self.SUM, prec)
        if op == "*":
            return self._product(args, prec, compact)
        if op == "/":
            num = self.render(args[0], self.TERM + 1, compact)
            den = self.render(args[1], self.TERM + 1, compact)
            joined = f"{num}/{den}" if compact or len(num) + len(den) <= 7 else f"{num} / {den}"
            return self._wrap(joined, self.TERM, prec)
        if op == "^":
            base = self.render(args[0], self.POW + 1, compact)
            exp_arg = args[1] if len(args) > 1 else None
            if (
                isinstance(exp_arg, (int, float))
                and not isinstance(exp_arg, bool)
                and float(exp_arg).is_integer()
            ):
                return base + sup(int(exp_arg))
            exponent = self.render(args[1], self.POW + 1, compact) if len(args) > 1 else "?"
            return self._wrap(f"{base}^{exponent}", self.POW, prec)
        if op == "index":
            base = self.render(args[0], self.ATOM, compact)
            idx = ", ".join(self.render(a, 0, compact=True) for a in args[1:])
            return f"{base}[{idx}]"
        if op == "apply_expression_template":
            name = node.get("name", "?")
            bindings = node.get("bindings") or {}
            if not bindings:
                return name
            rendered = ", ".join(self.render(v, 0, compact) for v in bindings.values())
            return f"{name}({rendered})"
        if op == "aggregate":
            expr, ranges, _ = self.render_aggregate(node)
            return f"[{expr} for {ranges}]"
        if op == "makearray":
            lines = []
            axes = [f"d{k + 1}" for k in range(len((node.get("regions") or [[]])[0]))]
            for label, role, value in self.render_makearray(node, axes):
                lines.append(f"{label} ({role}): {value}")
            return "makearray{ " + "; ".join(lines) + " }"
        if op == "D":
            return self._derivative(node, prec, compact)
        if op in ("<", ">", "<=", ">=", "==", "!="):
            pretty_op = {"<=": "≤", ">=": "≥", "==": "=", "!=": "≠"}.get(op, op)
            left = self.render(args[0], 31, compact)
            right = self.render(args[1], 31, compact)
            return self._wrap(f"{left} {pretty_op} {right}", 30, prec)

        # Unknown op: try the official display path for pure-math subtrees
        # (cos, exp, ...), else generic call syntax (recorded as a fallback).
        if self.official.available and official_safe(node):
            text = self.official.render(node)
            if text is not None:
                self.official_uses += 1
                return text
        self.structural_fallback_ops.add(str(op))
        extra = "".join(
            f", {k}: {node[k]}" for k in ("wrt", "dim") if k in node
        )
        rendered = ", ".join(self.render(a, 0, compact) for a in args)
        return f"{op}({rendered}{extra})"

    def render_makearray(self, node, axes: list[str]):
        """-> [(region label, role, rendered value), ...] in declaration order."""
        out = []
        regions = node.get("regions") or []
        values = node.get("values") or []
        for region, value in zip(regions, values):
            label, role = self.render_region(region, axes)
            if isinstance(value, dict) and value.get("op") == "aggregate":
                expr, ranges, _ = self.render_aggregate(value)
                out.append((label, role, f"{expr}   for {ranges}"))
            else:
                out.append((label, role, self.render_top(value)))
        return out

    # -- helpers --------------------------------------------------------------

    def _wrap(self, s: str, node_prec: int, prec: int) -> str:
        return f"({s})" if node_prec < prec else s

    def _number(self, value, prec: int) -> str:
        if isinstance(value, float) and value.is_integer():
            value = int(value)
        if isinstance(value, float):
            # Float display convention: at most 12 significant digits, matching
            # the official earthsci_toolkit.display path (display.py
            # _format_number: f"{num:.12g}" for regular-magnitude floats), so
            # official-rendered and structurally-rendered pages agree.
            # Presentation only — the underlying .esm bytes are untouched.
            if 1e-4 <= abs(value) < 1e5:
                text = f"{value:.12g}"
                if "." in text:
                    text = text.rstrip("0").rstrip(".")
            else:
                text = f"{value:.12g}"
        else:
            text = str(value)
        if text.startswith("-"):
            text = MINUS + text[1:]
            if prec > self.TERM:
                return f"({text})"
        return text

    def _negated(self, node):
        """Return a node equal to -node when the sign is syntactically evident."""
        if isinstance(node, (int, float)) and not isinstance(node, bool) and node < 0:
            return -node
        if isinstance(node, dict):
            if node.get("op") == "-" and len(node.get("args") or []) == 1:
                return node["args"][0]
            if node.get("op") == "*":
                args = node.get("args") or []
                if args and isinstance(args[0], (int, float)) and not isinstance(args[0], bool) and args[0] < 0:
                    first = -args[0]
                    rest = ([first] if first != 1 or len(args) == 1 else []) + list(args[1:])
                    return rest[0] if len(rest) == 1 else {"op": "*", "args": rest}
        return None

    def _sum(self, args, prec: int, compact: bool) -> str:
        if not args:
            return "0"
        text = self.render(args[0], self.SUM + 1, compact)
        for a in args[1:]:
            neg = self._negated(a)
            if neg is not None:
                sep = MINUS if compact else f" {MINUS} "
                text += sep + self.render(neg, self.SUM + 2, compact)
            else:
                text += ("+" if compact else " + ") + self.render(a, self.SUM + 1, compact)
        return self._wrap(text, self.SUM, prec)

    def _product(self, args, prec: int, compact: bool) -> str:
        if not args:
            return "1"
        if len(args) > 1 and all(a == args[0] for a in args[1:]):
            base = self.render(args[0], self.POW + 1, compact)
            return base + sup(len(args))
        prefix = ""
        if isinstance(args[0], (int, float)) and not isinstance(args[0], bool) and args[0] < 0:
            prefix = MINUS
            first = -args[0]
            args = ([first] if first != 1 or len(args) == 1 else []) + list(args[1:])
        # Every factor at TERM+1 so embedded divisions/sums parenthesize:
        # (a/b)·c, not the ambiguous a/b·c.
        parts = [self.render(a, self.TERM + 1, compact) for a in args]
        return self._wrap(prefix + DOT.join(parts), self.TERM, prec)

    def _derivative(self, node, prec: int, compact: bool) -> str:
        # Collect the chain of nested D nodes with a common wrt.
        order, inner, wrt = 0, node, node.get("wrt")
        while isinstance(inner, dict) and inner.get("op") == "D" and inner.get("wrt") == wrt:
            order += 1
            inner = (inner.get("args") or [None])[0]
        arg = self.render(inner, self.ATOM, compact)
        if order == 1:
            return f"{PARTIAL}{arg}/{PARTIAL}{wrt}"
        return f"{PARTIAL}{sup(order)}{arg}/{PARTIAL}{wrt}{sup(order)}"


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


def convergence_section(case: dict, renderer: MathRenderer, plots_root: Path) -> list[str]:
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
    stencil_templates: dict,
    renderer: MathRenderer,
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
    axes = axes_tag.split(",") if axes_tag else []
    if not axes:
        warn(f"{rel(rule_path)}: no axes: tag — using generic axis labels")

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
        pretty = renderer.render_top(match)
        priority_note = f" at priority {template['priority']}" if "priority" in template else ""
        out.append(f"**Rewrites** `{md_cell(pretty)}`{priority_note} — match pattern:")
        out.append("")
        out.append("```json")
        out.append(json.dumps(match, ensure_ascii=False))
        out.append("```")
        out.append("")
    else:
        warn(f"{rel(rule_path)}: rule template has no match pattern")

    # Discretization: the makearray region table (interior vs boundary faces).
    body = template.get("body")
    if isinstance(body, dict) and body.get("op") == "makearray":
        out.append(
            "**Discretization** — one `makearray` tiling the output axes; boundary "
            "conditions live in the face regions of this same rule (esm-spec §9.6.8; "
            "later regions overwrite earlier, §4.3.2):"
        )
        out.append("")
        out.append("| Region | Treatment | Value |")
        out.append("|---|---|---|")
        regions = body.get("regions") or []
        values = body.get("values") or []
        for region, value in zip(regions, values):
            label, role = renderer.render_region(region, axes)
            cell = render_region_value(value, stencil_templates, renderer)
            role_html = f'<span class="region-role {"face" if "face" in role else "interior"}">{role}</span>'
            out.append(f"| `{md_cell(label)}` | {role_html} | {md_cell(cell)} |")
        out.append("")
    elif body is not None:
        out.append("**Definition:**")
        out.append("")
        out.append(math_block(renderer.render_top(body)))
        out.append("")

    description = metadata.get("description")
    if description:
        out.append(description)
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
            out += convergence_section(case, renderer, plots_root)

    return out


def render_region_value(value, stencil_templates: dict, renderer: MathRenderer) -> str:
    """Pretty cell for a makearray region value; expands imported stencils."""
    if isinstance(value, dict) and value.get("op") == "apply_expression_template":
        name = value.get("name")
        applied = renderer.render(value, 0)
        stencil = stencil_templates.get(name)
        if stencil is not None:
            body = stencil.get("body")
            if isinstance(body, dict) and body.get("op") == "aggregate":
                expr, ranges, idx = renderer.render_aggregate(body)
                idx_s = ",".join(idx)
                return f"`{applied}[{idx_s}] = {expr}`   for `{ranges}`"
            return f"`{applied} = {renderer.render_top(body)}`"
        return f"`{applied}`"
    if isinstance(value, dict) and value.get("op") == "aggregate":
        expr, ranges, _ = renderer.render_aggregate(value)
        return f"`{expr}`   for `{ranges}`"
    return f"`{renderer.render_top(value)}`"


def stencil_section(stencil_path: Path, doc: dict, renderer: MathRenderer) -> list[str]:
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
    body = template.get("body")
    params = template.get("params") or []
    call = f"{name}({', '.join(params)})" if params else name
    if isinstance(body, dict) and body.get("op") == "aggregate":
        expr, ranges, idx = renderer.render_aggregate(body)
        idx_s = ",".join(idx)
        out.append(math_block(f"{call}[{idx_s}] = {expr}    for {ranges}"))
    elif body is not None:
        out.append(math_block(f"{call} = {renderer.render_top(body)}"))
    out.append("")
    description = metadata.get("description")
    if description:
        out.append(description)
        out.append("")
    out += references_md(metadata)
    return out


def grid_page(
    grid_dir: Path,
    renderer: MathRenderer,
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

    # Templates defined by this grid's stencils, for in-place expansion on rules.
    stencil_templates = {}
    for path, sdoc in stencil_docs:
        for tname, tmpl in (sdoc.get("expression_templates") or {}).items():
            stencil_templates[tname] = tmpl

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
        lines.append(metadata["description"])
        lines.append("")

    metaparameters = doc.get("metaparameters") or {}
    if metaparameters:
        lines += ["## Metaparameters", "", "| Name | Type | Default | Description |", "|---|---|---|---|"]
        for pname, spec in metaparameters.items():
            default = spec.get("default", "—")
            lines.append(
                f"| `{pname}` | {spec.get('type', '')} | `{default}` | {md_cell(spec.get('description', ''))} |"
            )
        lines += [
            "",
            "All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the "
            "import edge (or through the loader API) and the same files serve every resolution.",
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
        lines += ["## Geometry templates", "", "| Template | Definition | Description |", "|---|---|---|"]
        for tname, tmpl in templates.items():
            body = tmpl.get("body")
            if isinstance(body, dict) and body.get("op") == "aggregate":
                expr, ranges, idx = renderer.render_aggregate(body)
                idx_s = ",".join(idx)
                definition = f"`{tname}[{idx_s}] = {expr}` for `{ranges}`"
            else:
                definition = f"`{tname} = {renderer.render_top(body)}`"
            lines.append(f"| `{tname}` | {md_cell(definition)} | {md_cell(tmpl.get('description', ''))} |")
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
            lines += stencil_section(path, sdoc, renderer)

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
                path, rdoc, stencil_templates, renderer, conformance_cases, problem_cache, plots_root
            )
    else:
        lines += ["No rules published for this grid yet.", ""]

    return "\n".join(lines) + "\n"


def cross_grid_section_page(directory: Path, title: str, kind: str, renderer: MathRenderer) -> str:
    """Section page (_index.md) for regridding/ or reprojection/."""
    tag_prefix = {"regrid": "method:", "reproject": "crs:"}[kind]
    tag_label = {"regrid": "method", "reproject": "crs"}[kind]
    intro = {
        "regrid": "Conservative and interpolating remap expressions (coupling transforms). "
        "Each entry is a pure ESM template library — declarative expressions, no kernels.",
        "reproject": "Coordinate-transform template fragments (longitude-latitude, Lambert "
        "conformal, …) for use inside regridding and coupling expressions.",
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
    ]
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
            lines.append(metadata["description"])
            lines.append("")
        templates = doc.get("expression_templates") or {}
        if templates:
            lines.append(f"<details><summary>Templates ({len(templates)})</summary>")
            lines.append("")
            for tname, tmpl in templates.items():
                params = tmpl.get("params") or []
                call = f"{tname}({', '.join(params)})" if params else tname
                lines.append(f"**`{call}`**")
                if tmpl.get("description"):
                    lines.append("")
                    lines.append(tmpl["description"])
                lines.append("")
                body = tmpl.get("body")
                if isinstance(body, dict) and body.get("op") == "makearray":
                    axes = [f"d{k + 1}" for k in range(len((body.get("regions") or [[]])[0]))]
                    rows = renderer.render_makearray(body, axes)
                    text = "\n".join(f"{label}  ({role}):\n    {value}" for label, role, value in rows)
                    lines.append(math_block(text))
                elif isinstance(body, dict) and body.get("op") == "aggregate":
                    expr, ranges, idx = renderer.render_aggregate(body)
                    idx_s = ",".join(idx)
                    lines.append(math_block(f"{call}[{idx_s}] = {expr}    for {ranges}"))
                elif body is not None:
                    lines.append(math_block(f"{call} = {renderer.render_top(body)}"))
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
    official = OfficialDisplay()
    if official.available:
        print(f"[info] official ESS display path: {official.mode}")
    else:
        warn(
            "official ESS display path unavailable (earthsci_toolkit not importable and no "
            "ESS toolkit venv found) — all math rendered structurally"
        )
    renderer = MathRenderer(official)

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
            page = grid_page(grid_dir, renderer, conformance_cases, problem_cache, plots_root)
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
            cross_grid_section_page(directory, title, kind, renderer), encoding="utf-8"
        )

    official.close()

    print(f"[info] generated {n_grids} grid page(s) into {grids_out}")
    if renderer.official_uses:
        print(f"[info] official display path rendered {renderer.official_uses} expression(s)")
    if renderer.structural_fallback_ops:
        print(
            "[info] ops rendered with generic structural fallback (no official rendering): "
            + ", ".join(sorted(renderer.structural_fallback_ops))
        )


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
