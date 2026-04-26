# EarthSciDiscretizations Hugo doc site

This directory holds the Hugo site that mirrors EarthSciModels' component
catalog, but for ESD's discretization rules and grid families.

## Layout

```
docs/
├── hugo.toml                    site config (theme-less; layouts inline)
├── content/
│   ├── _index.md                home page
│   ├── grids/
│   │   ├── _index.md            section index
│   │   └── <family>.md          one page per grid family
│   └── rules/
│       ├── _index.md            section index
│       └── <rule>.md            one page per seeded rule
├── layouts/
│   ├── _default/{baseof,list,single,terms}.html
│   └── partials/{head,header,footer}.html
└── static/
    ├── css/site.css
    └── plots/
        ├── grids/<family>.png
        └── rules/<rule>-{stencil,convergence}.png
```

The Documenter.jl site (under `src/` with `make.jl`) is a sibling — Hugo
reads only `content/`, `layouts/`, `static/`, and `data/`, so the two
docsets coexist without conflict.

## Build

```bash
# 1. Render plot artifacts into static/plots/.
python3 ../tools/render_doc_plots.py

# 2. Build the Hugo site.
hugo --source . --minify --destination public
```

Open `public/index.html` for the result. CI must run **step 1 before
step 2** so the static plot files are in place when Hugo copies
`static/`.

## Plot rendering

`tools/render_doc_plots.py` is a self-contained matplotlib script — no
cartopy, no Julia. It produces:

- 7 grid-family visualizations (cartesian, latlon, cubed_sphere, mpas,
  duo, vertical, arakawa)
- 11 rule stencil / coefficient diagrams
- 3 convergence plots — for the rules whose Layer-B fixtures are
  currently producing (`centered_2nd_uniform`,
  `centered_2nd_uniform_vertical`, `upwind_1st`).

Convergence plots for the other 8 rules are **gated on in-flight ESS
harness extensions** in
`EarthSciSerialization/packages/EarthSciSerialization.jl/src/mms_evaluator.jl`
(2D dispatch + per-cell metric callables, MPAS unstructured support,
sub-stencil targeting + parabola pass for PPM, nonlinear reconstruction
for WENO5, stagger-aware sampling for Arakawa-C divergence, monotonicity
/ TVD harness for the limiters). Each affected page renders a *pending*
callout in its convergence section until the corresponding harness
extension lands and the fixture's `applicable` flag flips to true; at
that point `render_doc_plots.py` will be extended to consume the walker
output directly.

## Adding a new rule or grid family

1. Drop the rule JSON / grid trait in the usual place
   (`discretizations/<family>/<name>.json`,
   `src/grids/<family>.jl`).
2. Add a stencil renderer (and convergence renderer if applicable) to
   `tools/render_doc_plots.py`.
3. Add a markdown page under `content/rules/` or `content/grids/`.
4. Re-run the two build steps.
