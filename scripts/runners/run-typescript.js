#!/usr/bin/env node
// run-typescript.js — the official TypeScript conformance runner for
// EarthSciDiscretizations.
//
// THIN WRAPPER (AGENTS.md §2, single pathway): every byte this script emits
// comes from official earthsci-ast-ts entry points driven over the one
// canonical pipeline (.esm → parse → §9.7 import/metaparameter resolution →
// §9.6.3 rewrite fixpoint). No evaluator, rule engine, or numeric kernel
// lives here — only argument marshalling and JSON I/O.
//
// CLI (EarthSciAST CONFORMANCE_SPEC.md §4.2):
//   node run-typescript.js --output-dir <path> [--categories ast]
//        [--files <manifest.json>[,…]] [--verbose]
//
// Scope: the TypeScript binding is the rewrite half of the standard (full
// §9.6.3/§9.7 lowering, no makearray simulator — see the ESD manifests'
// scope_excluded notes), so this runner implements the `ast` category only.
//
// Outputs (§4.4 layout, mirroring scripts/runners/run-julia.jl):
//   <output-dir>/typescript_ast_results.json
//   <output-dir>/typescript_summary.json
//   <output-dir>/ast/<case>/expanded.golden.json   canonical post-lowering
//                                                  bytes (byte-compared vs
//                                                  the golden)
//
// Environment: requires the built earthsci-ast-ts dist
// ($ESS_ROOT/pkg/earthsci-ast-ts → npm ci && npm run build).

'use strict'

const fs = require('fs')
const path = require('path')

const ESD_ROOT = path.resolve(__dirname, '..', '..')
const ESS_ROOT = process.env.ESS_ROOT
  ? path.resolve(process.env.ESS_ROOT)
  : path.resolve(ESD_ROOT, '..', 'EarthSciAST')
if (!fs.existsSync(path.join(ESS_ROOT, 'esm-schema.json'))) {
  console.error(`error: EarthSciAST not found at '${ESS_ROOT}'; ` +
    'set ESS_ROOT or clone it as a sibling checkout (scripts/ess-locate.sh contract)')
  process.exit(1)
}
const DIST = path.join(ESS_ROOT, 'pkg', 'earthsci-ast-ts', 'dist', 'cjs', 'index.js')
if (!fs.existsSync(DIST)) {
  console.error(`error: earthsci-ast-ts dist not built at '${DIST}' ` +
    '(cd pkg/earthsci-ast-ts && npm ci && npm run build)')
  process.exit(1)
}
const { resolveTemplateMachinery, lowerExpressionTemplates } = require(DIST)

const CONF = path.join(ESD_ROOT, 'tests', 'conformance')

// ---------------------------------------------------------------------------
// Canonical golden writer — byte-for-byte the Julia reference runner's scheme
// (scripts/runners/run-julia.jl): object keys sorted (code-unit order), arrays
// in order, 2-space indent, trailing newline, scalars as JSON3 renders them
// after JSON3.read's numeric normalization (integral floats → integers;
// non-integral floats via Ryu shortest with Julia's positional/scientific
// switch at 1e-5 / 1e6).
// ---------------------------------------------------------------------------

function juliaFloatRepr (v) {
  if (!Number.isFinite(v)) throw new Error(`non-finite float ${v} is not valid JSON`)
  let s = String(Math.abs(v)) // shortest round-trip digits
  const neg = v < 0
  let mant = s
  let e = 0
  const em = s.match(/^([^eE]+)[eE]([+-]?\d+)$/)
  if (em) { mant = em[1]; e = parseInt(em[2], 10) }
  const dot = mant.indexOf('.')
  const ip = dot === -1 ? mant : mant.slice(0, dot)
  const fp = dot === -1 ? '' : mant.slice(dot + 1)
  let digits = (ip + fp).replace(/^0+/, '')
  let e10
  if (ip !== '' && ip !== '0' && !/^0+$/.test(ip)) {
    e10 = e + ip.length - 1
  } else {
    const lead = fp.length - fp.replace(/^0+/, '').length
    e10 = e - lead - 1
  }
  digits = digits.replace(/0+$/, '') || '0'
  const sign = neg ? '-' : ''
  if (e10 >= -4 && e10 <= 5) {
    if (e10 < 0) return sign + '0.' + '0'.repeat(-e10 - 1) + digits
    if (digits.length <= e10 + 1) {
      return sign + digits + '0'.repeat(e10 + 1 - digits.length) + '.0'
    }
    return sign + digits.slice(0, e10 + 1) + '.' + digits.slice(e10 + 1)
  }
  const frac = digits.slice(1) || '0'
  return sign + digits[0] + '.' + frac + 'e' + String(e10)
}

const INT64_MAX = 9223372036854775807n

function scalar (x) {
  if (x === true) return 'true'
  if (x === false) return 'false'
  if (x === null) return 'null'
  if (typeof x === 'number') {
    // JSON3.read numeric normalization: integral values parse as Int64.
    if (Number.isInteger(x) && Math.abs(x) <= Number.MAX_SAFE_INTEGER) return String(x)
    if (Number.isInteger(x) && BigInt(x) <= INT64_MAX && BigInt(x) >= -INT64_MAX - 1n) {
      return BigInt(x).toString()
    }
    return juliaFloatRepr(x)
  }
  if (typeof x === 'string') return JSON.stringify(x)
  throw new Error(`unsupported scalar ${typeof x}`)
}

function writeSorted (out, x, indent) {
  const pad = '  '.repeat(indent)
  const pad1 = '  '.repeat(indent + 1)
  if (Array.isArray(x)) {
    if (x.length === 0) { out.push('[]'); return }
    out.push('[\n')
    x.forEach((v, i) => {
      out.push(pad1)
      writeSorted(out, v, indent + 1)
      out.push(i < x.length - 1 ? ',\n' : '\n')
    })
    out.push(pad + ']')
  } else if (x !== null && typeof x === 'object') {
    const ks = Object.keys(x).sort()
    if (ks.length === 0) { out.push('{}'); return }
    out.push('{\n')
    ks.forEach((k, i) => {
      out.push(pad1 + JSON.stringify(k) + ': ')
      writeSorted(out, x[k], indent + 1)
      out.push(i < ks.length - 1 ? ',\n' : '\n')
    })
    out.push(pad + '}')
  } else {
    out.push(scalar(x))
  }
}

function canonicalBytes (doc) {
  const out = []
  writeSorted(out, doc, 0)
  out.push('\n')
  return Buffer.from(out.join(''), 'utf8')
}

// ---------------------------------------------------------------------------
// Manifest discovery + the post-lowering structural gates (mirrors
// run-julia.jl's unlowered_violations).
// ---------------------------------------------------------------------------

function discoverManifests (category, files) {
  if (files.length > 0) {
    return files
      .map(f => path.isAbsolute(f) ? f : path.join(ESD_ROOT, f))
      .map(p => [path.dirname(path.resolve(p)), JSON.parse(fs.readFileSync(p, 'utf8'))])
      .filter(([, m]) => m.category === category)
  }
  const dir = path.join(CONF, category)
  if (!fs.existsSync(dir)) return []
  return fs.readdirSync(dir).sort()
    .map(c => path.join(dir, c, 'manifest.json'))
    .filter(mp => fs.existsSync(mp))
    .map(mp => [path.dirname(mp), JSON.parse(fs.readFileSync(mp, 'utf8'))])
}

// The prose-free structural gates (mirrors run-julia.jl).
//
// NOTE (ESS 0.9.0, esd:duo migration): surviving `apply_expression_template`
// nodes are NO LONGER a violation. Since the Option B reference-preserving
// switch (EarthSciAST commit aff96f29, esm 0.9.0) the loader does not inline
// match-less template bodies — a reference "denotes its expansion" and the
// engine treats it as a leaf (esm-spec §9.6.4). Because `match` patterns may
// not contain `apply_expression_template` (§9.6.1), every surviving apply is a
// resolved match-less leaf, not an un-fired rule; a dangling apply to an
// unknown template is already rejected upstream by resolveTemplateMachinery.
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
function unloweredViolations (node, acc) {
  if (Array.isArray(node)) {
    node.forEach(v => unloweredViolations(v, acc))
  } else if (node !== null && typeof node === 'object') {
    const op = String(node.op || '')
    if (op === 'grad' || op === 'div' || op === 'laplacian') acc.push(`unlowered ${op}`)
    if (op === 'D' && String(node.wrt || 't') !== 't') {
      acc.push(`unlowered spatial D (wrt=${node.wrt})`)
    }
    for (const [k, v] of Object.entries(node)) {
      // `metadata`: prose may cite op names.
      // `expression_templates`: see the note above this function.
      if (k === 'metadata' || k === 'expression_templates') continue
      unloweredViolations(v, acc)
    }
  }
  return acc
}

function runAst (outputDir, files, verbose) {
  const cases = {}
  for (const [caseDir, manifest] of discoverManifests('ast', files)) {
    const kase = manifest.case
    const rec = { case: kase, status: 'ok' }
    try {
      const fixture = path.join(caseDir, manifest.fixture)
      const raw = JSON.parse(fs.readFileSync(fixture, 'utf8'))
      const resolved = resolveTemplateMachinery(raw, path.dirname(fixture))
      const doc = lowerExpressionTemplates(resolved !== null ? resolved : raw)
      const violations = []
      for (const model of Object.values(doc.models || {})) {
        const stripped = Object.fromEntries(
          Object.entries(model).filter(
            ([k]) => k !== 'metadata' && k !== 'expression_templates'))
        unloweredViolations(stripped, violations)
      }
      if (violations.length > 0) {
        throw new Error('post-lowering gate failed: ' +
          [...new Set(violations)].sort().join(', '))
      }
      const bytes = canonicalBytes(doc)
      const artDir = path.join(outputDir, 'ast', kase)
      fs.mkdirSync(artDir, { recursive: true })
      const art = path.join(artDir, 'expanded.golden.json')
      fs.writeFileSync(art, bytes)
      rec.artifact = path.relative(outputDir, art)
      rec.bytes = bytes.length
      rec.golden = path.relative(ESD_ROOT, path.join(caseDir, manifest.golden))
    } catch (err) {
      rec.status = 'error'
      rec.message = String(err && err.message ? err.message : err)
    }
    if (verbose) console.log(`  ast/${kase}: ${rec.status}`)
    cases[kase] = rec
  }
  return { binding: 'typescript', category: 'ast', cases }
}

// ---------------------------------------------------------------------------
// Main (§4.2 CLI).
// ---------------------------------------------------------------------------

function main () {
  const argv = process.argv.slice(2)
  let outputDir = null
  let categories = []
  let files = []
  let verbose = false
  for (let i = 0; i < argv.length; i++) {
    switch (argv[i]) {
      case '--output-dir': outputDir = argv[++i]; break
      case '--categories': categories = argv[++i].split(',').filter(Boolean); break
      case '--files': files = argv[++i].split(',').filter(Boolean); break
      case '--verbose': verbose = true; break
      default:
        console.error(`error: unknown argument '${argv[i]}'`)
        return 2
    }
  }
  if (!outputDir) {
    console.error('error: --output-dir is required (CONFORMANCE_SPEC.md §4.2)')
    return 2
  }
  if (categories.length === 0) categories = ['ast']
  if (files.length > 0) {
    const fileCats = new Set(files.map(f =>
      JSON.parse(fs.readFileSync(path.isAbsolute(f) ? f : path.join(ESD_ROOT, f),
        'utf8')).category))
    categories = categories.filter(c => fileCats.has(c))
  }
  fs.mkdirSync(outputDir, { recursive: true })
  const summary = {
    binding: 'typescript',
    runner: 'scripts/runners/run-typescript.js',
    categories: {}
  }
  let exitCode = 0
  for (const cat of categories) {
    if (cat !== 'ast') {
      console.log(`typescript/${cat}: unsupported category (skipping; ` +
        'rewrite-only port — ast only)')
      continue
    }
    if (verbose) console.log(`category: ${cat}`)
    const res = runAst(outputDir, files, verbose)
    fs.writeFileSync(path.join(outputDir, `typescript_${cat}_results.json`),
      canonicalBytes(res))
    const recs = Object.values(res.cases)
    const nBad = recs.filter(r => r.status !== 'ok').length
    if (nBad > 0) exitCode = 1
    summary.categories[cat] = { cases: recs.length, failures: nBad }
    console.log(`typescript/${cat}: ${recs.length} case(s), ${nBad} failure(s)`)
  }
  fs.writeFileSync(path.join(outputDir, 'typescript_summary.json'),
    canonicalBytes(summary))
  return exitCode
}

process.exit(main())
