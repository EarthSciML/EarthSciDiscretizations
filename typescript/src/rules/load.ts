/**
 * Loaders for ESS-emitted rule JSON files in the discretization catalog.
 */

import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, basename, extname } from "node:path";

import type { Rule } from "./types.js";

interface RawRuleFile {
  discretizations: Record<string, Record<string, unknown>>;
}

/**
 * Parse a single rule JSON file. The on-disk layout is:
 *
 *   { "discretizations": { "<rule_name>": { applies_to, grid_family, stencil, ... } } }
 *
 * The family is taken from the immediate parent directory name (e.g.,
 * `discretizations/finite_difference/centered_2nd_uniform_latlon.json`
 * → family `finite_difference`).
 */
export function parseRuleFile(path: string, content: string): Rule[] {
  const parsed = JSON.parse(content) as RawRuleFile;
  if (!parsed.discretizations) {
    throw new Error(`rule file ${path} missing 'discretizations' top-level key`);
  }
  const family = basename(path.replace(/[/\\][^/\\]+$/, ""));
  const out: Rule[] = [];
  for (const [name, body] of Object.entries(parsed.discretizations)) {
    out.push({ name, family, ...body } as unknown as Rule);
  }
  return out;
}

/** Load and parse one rule JSON file. */
export function loadRuleFile(path: string): Rule[] {
  return parseRuleFile(path, readFileSync(path, "utf8"));
}

/**
 * Walk a discretization catalog (e.g., `<repo>/discretizations`) and collect
 * every `<family>/<rule>.json` entry. Sibling directories named after the rule
 * (which hold fixtures) are ignored.
 */
export function loadCatalog(catalogDir: string): Rule[] {
  const rules: Rule[] = [];
  for (const familyEntry of readdirSync(catalogDir)) {
    const familyPath = join(catalogDir, familyEntry);
    if (!statSync(familyPath).isDirectory()) continue;
    for (const child of readdirSync(familyPath)) {
      if (extname(child) !== ".json") continue;
      const childPath = join(familyPath, child);
      if (!statSync(childPath).isFile()) continue;
      rules.push(...loadRuleFile(childPath));
    }
  }
  return rules;
}

/** Find a rule by name across the loaded catalog. Returns `undefined` if absent. */
export function findRule(rules: Rule[], name: string): Rule | undefined {
  return rules.find((r) => r.name === name);
}
