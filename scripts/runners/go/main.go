// run-go — the official Go conformance runner for EarthSciDiscretizations
// (invoked through scripts/runners/run-go.sh, which resolves the ESS module).
//
// THIN WRAPPER (AGENTS.md §2, single pathway): every byte this program emits
// comes from the official earthsci-ast-go entry point
// esm.ResolveAndLowerReferencePreserving — the raw §9.7 pipeline (.esm → parse
// → import/metaparameter resolution → §9.6.3 rewrite fixpoint), stopping at the
// Option-B image. No evaluator, rule engine, or numeric kernel lives here —
// only argument marshalling and JSON I/O.
//
// It is deliberately NOT esm.ResolveAndLower, which appends Expand-at-build.
// The other four runners all stop at the Option-B image, so expanding here made
// Go's artifact structurally incomparable with theirs; and expansion multiplies
// out a reference-dense document — every duo case OOM-killed the runner past
// 5 GB, so Go produced no artifact for that family at all.
//
// CLI (EarthSciAST CONFORMANCE_SPEC.md §4.2):
//
//	run-go --output-dir <path> [--categories ast]
//	       [--files <manifest.json>[,…]] [--verbose]
//
// Scope: the Go binding is the rewrite half of the standard (full §9.6.3/§9.7
// lowering, no makearray simulator — see the ESD manifests' scope_excluded
// notes), so this runner implements the `ast` category only.
//
// Outputs (§4.4 layout, mirroring scripts/runners/run-julia.jl):
//
//	<output-dir>/go_ast_results.json
//	<output-dir>/go_summary.json
//	<output-dir>/ast/<case>/expanded.golden.json  canonical post-lowering bytes
//
// The environment variable ESD_ROOT points at the EarthSciDiscretizations
// checkout (set by run-go.sh).
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	esm "github.com/EarthSciML/EarthSciAST/pkg/earthsci-ast-go/pkg/esm"
)

// ---------------------------------------------------------------------------
// Canonical golden writer — byte-for-byte the Julia reference runner's scheme
// (scripts/runners/run-julia.jl): object keys sorted (byte order), arrays in
// order, 2-space indent, trailing newline, scalars as JSON3 renders them after
// JSON3.read's numeric normalization (integral floats → Int64; non-integral
// floats via Ryu shortest with Julia's positional/scientific switch at
// 1e-5 / 1e6; strings UTF-8 with the JSON3 escape table).
// ---------------------------------------------------------------------------

// juliaFloatRepr renders a non-integral finite float exactly as Julia's
// JSON3.write (Base.Ryu shortest digits; scientific iff the decimal exponent
// is < -4 or >= 6; scientific mantissa always carries a fraction digit;
// exponent unpadded, no '+').
func juliaFloatRepr(v float64) (string, error) {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return "", fmt.Errorf("non-finite float %v is not valid JSON", v)
	}
	s := strconv.FormatFloat(math.Abs(v), 'e', -1, 64) // d[.ddd]e±XX shortest
	mantExp := strings.SplitN(s, "e", 2)
	mant := mantExp[0]
	e, err := strconv.Atoi(mantExp[1])
	if err != nil {
		return "", err
	}
	digits := strings.Replace(mant, ".", "", 1)
	digits = strings.TrimRight(digits, "0")
	if digits == "" {
		digits = "0"
	}
	// FormatFloat 'e' mantissa is d[.ddd] with one leading digit, so the
	// decimal exponent of the first significant digit is e itself.
	e10 := e
	sign := ""
	if v < 0 {
		sign = "-"
	}
	if e10 >= -4 && e10 <= 5 { // positional
		if e10 < 0 {
			return sign + "0." + strings.Repeat("0", -e10-1) + digits, nil
		}
		if len(digits) <= e10+1 {
			return sign + digits + strings.Repeat("0", e10+1-len(digits)) + ".0", nil
		}
		return sign + digits[:e10+1] + "." + digits[e10+1:], nil
	}
	frac := "0"
	if len(digits) > 1 {
		frac = digits[1:]
	}
	return sign + digits[:1] + "." + frac + "e" + strconv.Itoa(e10), nil
}

// escapeString writes a JSON string with the JSON3 escape table: short
// escapes for \" \\ \b \f \n \r \t, \u00xx for other control characters,
// everything else (including non-ASCII) as raw UTF-8.
func escapeString(s string) string {
	var b strings.Builder
	b.WriteByte('"')
	for _, r := range s {
		switch r {
		case '"':
			b.WriteString(`\"`)
		case '\\':
			b.WriteString(`\\`)
		case '\b':
			b.WriteString(`\b`)
		case '\f':
			b.WriteString(`\f`)
		case '\n':
			b.WriteString(`\n`)
		case '\r':
			b.WriteString(`\r`)
		case '\t':
			b.WriteString(`\t`)
		default:
			if r < 0x20 {
				fmt.Fprintf(&b, `\u%04x`, r)
			} else {
				b.WriteRune(r)
			}
		}
	}
	b.WriteByte('"')
	return b.String()
}

func scalar(x interface{}) (string, error) {
	switch v := x.(type) {
	case nil:
		return "null", nil
	case bool:
		if v {
			return "true", nil
		}
		return "false", nil
	case string:
		return escapeString(v), nil
	case int64:
		return strconv.FormatInt(v, 10), nil
	case float64:
		return juliaFloatRepr(v)
	case json.Number:
		return "", fmt.Errorf("json.Number %q must be normalized before writing", v)
	}
	return "", fmt.Errorf("unsupported scalar %T", x)
}

// norm mirrors JSON3.read's numeric normalization: every json.Number that is
// integral and fits in Int64 becomes int64 (0.0 → 0), anything else float64.
func norm(x interface{}) (interface{}, error) {
	switch v := x.(type) {
	case map[string]interface{}:
		out := make(map[string]interface{}, len(v))
		for k, c := range v {
			n, err := norm(c)
			if err != nil {
				return nil, err
			}
			out[k] = n
		}
		return out, nil
	case []interface{}:
		out := make([]interface{}, len(v))
		for i, c := range v {
			n, err := norm(c)
			if err != nil {
				return nil, err
			}
			out[i] = n
		}
		return out, nil
	case json.Number:
		if i, err := v.Int64(); err == nil {
			return i, nil
		}
		f, err := v.Float64()
		if err != nil {
			return nil, fmt.Errorf("unparsable number %q", v)
		}
		if f == math.Trunc(f) && f >= math.MinInt64 && f < math.MaxInt64 {
			return int64(f), nil
		}
		return f, nil
	}
	return x, nil
}

func writeSorted(b *strings.Builder, x interface{}, indent int) error {
	pad := strings.Repeat("  ", indent)
	pad1 := strings.Repeat("  ", indent+1)
	switch v := x.(type) {
	case map[string]interface{}:
		if len(v) == 0 {
			b.WriteString("{}")
			return nil
		}
		ks := make([]string, 0, len(v))
		for k := range v {
			ks = append(ks, k)
		}
		sort.Strings(ks)
		b.WriteString("{\n")
		for i, k := range ks {
			b.WriteString(pad1 + escapeString(k) + ": ")
			if err := writeSorted(b, v[k], indent+1); err != nil {
				return err
			}
			if i < len(ks)-1 {
				b.WriteString(",")
			}
			b.WriteString("\n")
		}
		b.WriteString(pad + "}")
	case []interface{}:
		if len(v) == 0 {
			b.WriteString("[]")
			return nil
		}
		b.WriteString("[\n")
		for i, c := range v {
			b.WriteString(pad1)
			if err := writeSorted(b, c, indent+1); err != nil {
				return err
			}
			if i < len(v)-1 {
				b.WriteString(",")
			}
			b.WriteString("\n")
		}
		b.WriteString(pad + "]")
	default:
		s, err := scalar(x)
		if err != nil {
			return err
		}
		b.WriteString(s)
	}
	return nil
}

func canonicalBytes(doc interface{}) ([]byte, error) {
	n, err := norm(doc)
	if err != nil {
		return nil, err
	}
	var b strings.Builder
	if err := writeSorted(&b, n, 0); err != nil {
		return nil, err
	}
	b.WriteString("\n")
	return []byte(b.String()), nil
}

// ---------------------------------------------------------------------------
// Manifest discovery + the post-lowering structural gates (mirrors
// run-julia.jl's unlowered_violations).
// ---------------------------------------------------------------------------

type manifestEntry struct {
	dir      string
	manifest map[string]interface{}
}

func discoverManifests(esdRoot, category string, files []string) ([]manifestEntry, error) {
	var out []manifestEntry
	if len(files) > 0 {
		for _, f := range files {
			p := f
			if !filepath.IsAbs(p) {
				p = filepath.Join(esdRoot, p)
			}
			m, err := readJSONMap(p)
			if err != nil {
				return nil, err
			}
			if s, _ := m["category"].(string); s == category {
				abs, _ := filepath.Abs(p)
				out = append(out, manifestEntry{filepath.Dir(abs), m})
			}
		}
		return out, nil
	}
	dir := filepath.Join(esdRoot, "tests", "conformance", category)
	entries, err := os.ReadDir(dir)
	if err != nil {
		if os.IsNotExist(err) {
			return out, nil
		}
		return nil, err
	}
	names := make([]string, 0, len(entries))
	for _, e := range entries {
		names = append(names, e.Name())
	}
	sort.Strings(names)
	for _, name := range names {
		mp := filepath.Join(dir, name, "manifest.json")
		if _, err := os.Stat(mp); err != nil {
			continue
		}
		m, err := readJSONMap(mp)
		if err != nil {
			return nil, err
		}
		out = append(out, manifestEntry{filepath.Join(dir, name), m})
	}
	return out, nil
}

func readJSONMap(path string) (map[string]interface{}, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var m map[string]interface{}
	if err := json.Unmarshal(data, &m); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return m, nil
}

// unloweredViolations is the prose-free structural gate (mirrors run-julia.jl).
//
// NOTE (ESS 0.9.0, esd:duo migration): surviving apply_expression_template
// nodes are NO LONGER a violation. Since the Option B reference-preserving
// switch (EarthSciAST commit aff96f29, esm 0.9.0) the loader does not inline
// match-less template bodies — a reference "denotes its expansion" and the
// engine treats it as a leaf (esm-spec §9.6.4). Because match patterns may not
// contain apply_expression_template (§9.6.1), every surviving apply is a
// resolved match-less leaf, not an un-fired rule; a dangling apply to an
// unknown template is already rejected upstream by
// esm.ResolveAndLowerReferencePreserving (stopping before Expand skips no
// validation — the resolver, the fixpoint and the §9.6.9 call-site check all
// still reject).
// A model's retained `expression_templates` registry is NOT walked. Option
// B materializes each component's referenced templates into it, and a
// rule's `match` pattern is a PATTERN — `D(f, x)` there is what the rule
// fires ON, not an operator that survived lowering. §9.6.3 constraint 6
// scopes `unlowered_operator` to expressions reaching EVALUATION or
// COMPILATION (esm-spec: "fails at evaluation/compilation"); a match
// pattern reaches neither. Sweeping all 139 lowered ast fixtures: every
// rewrite-target op lives in a retained registry, ZERO in an evaluation
// position — so this clause has never once caught a real unlowered operator
// here, while failing every case that uses a spatial-D rule. Upstream
// already skips registries the same way (`_EXPR_TEMPLATES_SKIP` in Python,
// the `expression_templates` continue in Julia), and the real gate in
// tree_walk/compile.jl never walks one.
func unloweredViolations(node interface{}, acc *[]string) {
	switch v := node.(type) {
	case map[string]interface{}:
		op, _ := v["op"].(string)
		if op == "grad" || op == "div" || op == "laplacian" {
			*acc = append(*acc, "unlowered "+op)
		}
		if op == "D" {
			if wrt, ok := v["wrt"].(string); ok && wrt != "t" {
				*acc = append(*acc, "unlowered spatial D (wrt="+wrt+")")
			}
		}
		for k, c := range v {
			// "metadata": prose may cite op names.
			// "expression_templates": see the note above this function.
			if k == "metadata" || k == "expression_templates" {
				continue
			}
			unloweredViolations(c, acc)
		}
	case []interface{}:
		for _, c := range v {
			unloweredViolations(c, acc)
		}
	}
}

func uniqueSorted(in []string) []string {
	seen := map[string]bool{}
	var out []string
	for _, s := range in {
		if !seen[s] {
			seen[s] = true
			out = append(out, s)
		}
	}
	sort.Strings(out)
	return out
}

// ---------------------------------------------------------------------------
// ast — raw §9.7 expansion to canonical bytes + post-lowering gates.
// ---------------------------------------------------------------------------

func runAst(esdRoot, outputDir string, files []string, verbose bool) (map[string]interface{}, error) {
	cases := map[string]interface{}{}
	mans, err := discoverManifests(esdRoot, "ast", files)
	if err != nil {
		return nil, err
	}
	for _, me := range mans {
		caseName, _ := me.manifest["case"].(string)
		rec := map[string]interface{}{"case": caseName, "status": "ok"}
		err := func() error {
			fixture := filepath.Join(me.dir, me.manifest["fixture"].(string))
			data, err := os.ReadFile(fixture)
			if err != nil {
				return err
			}
			lowered, err := esm.ResolveAndLowerReferencePreserving(
				string(data), filepath.Dir(fixture), nil)
			if err != nil {
				return err
			}
			dec := json.NewDecoder(strings.NewReader(lowered))
			dec.UseNumber()
			var doc map[string]interface{}
			if err := dec.Decode(&doc); err != nil {
				return err
			}
			var violations []string
			if models, ok := doc["models"].(map[string]interface{}); ok {
				for _, model := range models {
					if mm, ok := model.(map[string]interface{}); ok {
						for k, v := range mm {
							if k == "metadata" || k == "expression_templates" {
								continue
							}
							unloweredViolations(v, &violations)
						}
					}
				}
			}
			if len(violations) > 0 {
				return fmt.Errorf("post-lowering gate failed: %s",
					strings.Join(uniqueSorted(violations), ", "))
			}
			bytes, err := canonicalBytes(doc)
			if err != nil {
				return err
			}
			artDir := filepath.Join(outputDir, "ast", caseName)
			if err := os.MkdirAll(artDir, 0o755); err != nil {
				return err
			}
			art := filepath.Join(artDir, "expanded.golden.json")
			if err := os.WriteFile(art, bytes, 0o644); err != nil {
				return err
			}
			relArt, _ := filepath.Rel(outputDir, art)
			relGold, _ := filepath.Rel(esdRoot,
				filepath.Join(me.dir, me.manifest["golden"].(string)))
			rec["artifact"] = relArt
			rec["bytes"] = int64(len(bytes))
			rec["golden"] = relGold
			return nil
		}()
		if err != nil {
			rec["status"] = "error"
			rec["message"] = err.Error()
		}
		if verbose {
			fmt.Printf("  ast/%s: %v\n", caseName, rec["status"])
		}
		cases[caseName] = rec
	}
	return map[string]interface{}{
		"binding": "go", "category": "ast", "cases": cases,
	}, nil
}

// ---------------------------------------------------------------------------
// Main (§4.2 CLI).
// ---------------------------------------------------------------------------

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(2)
	}
}

func run() error {
	esdRoot := os.Getenv("ESD_ROOT")
	if esdRoot == "" {
		return fmt.Errorf("ESD_ROOT is not set (run through scripts/runners/run-go.sh)")
	}
	args := os.Args[1:]
	outputDir := ""
	var categories, files []string
	verbose := false
	for i := 0; i < len(args); i++ {
		switch args[i] {
		case "--output-dir":
			i++
			outputDir = args[i]
		case "--categories":
			i++
			categories = splitCSV(args[i])
		case "--files":
			i++
			files = splitCSV(args[i])
		case "--verbose":
			verbose = true
		default:
			return fmt.Errorf("unknown argument %q", args[i])
		}
	}
	if outputDir == "" {
		return fmt.Errorf("--output-dir is required (CONFORMANCE_SPEC.md §4.2)")
	}
	if len(categories) == 0 {
		categories = []string{"ast"}
	}
	if len(files) > 0 {
		fileCats := map[string]bool{}
		for _, f := range files {
			p := f
			if !filepath.IsAbs(p) {
				p = filepath.Join(esdRoot, p)
			}
			m, err := readJSONMap(p)
			if err != nil {
				return err
			}
			if s, _ := m["category"].(string); s != "" {
				fileCats[s] = true
			}
		}
		var kept []string
		for _, c := range categories {
			if fileCats[c] {
				kept = append(kept, c)
			}
		}
		categories = kept
	}
	if err := os.MkdirAll(outputDir, 0o755); err != nil {
		return err
	}
	summary := map[string]interface{}{
		"binding":    "go",
		"runner":     "scripts/runners/run-go.sh",
		"categories": map[string]interface{}{},
	}
	exitCode := 0
	for _, cat := range categories {
		if cat != "ast" {
			fmt.Printf("go/%s: unsupported category (skipping; rewrite-only port — ast only)\n", cat)
			continue
		}
		if verbose {
			fmt.Printf("category: %s\n", cat)
		}
		res, err := runAst(esdRoot, outputDir, files, verbose)
		if err != nil {
			return err
		}
		bytes, err := canonicalBytes(res)
		if err != nil {
			return err
		}
		if err := os.WriteFile(filepath.Join(outputDir, "go_"+cat+"_results.json"),
			bytes, 0o644); err != nil {
			return err
		}
		cases := res["cases"].(map[string]interface{})
		nBad := 0
		for _, r := range cases {
			if r.(map[string]interface{})["status"] != "ok" {
				nBad++
			}
		}
		if nBad > 0 {
			exitCode = 1
		}
		summary["categories"].(map[string]interface{})[cat] = map[string]interface{}{
			"cases": int64(len(cases)), "failures": int64(nBad),
		}
		fmt.Printf("go/%s: %d case(s), %d failure(s)\n", cat, len(cases), nBad)
	}
	sumBytes, err := canonicalBytes(summary)
	if err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(outputDir, "go_summary.json"),
		sumBytes, 0o644); err != nil {
		return err
	}
	os.Exit(exitCode)
	return nil
}

func splitCSV(s string) []string {
	var out []string
	for _, p := range strings.Split(s, ",") {
		if p != "" {
			out = append(out, p)
		}
	}
	return out
}
