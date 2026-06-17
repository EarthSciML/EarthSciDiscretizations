using Documenter
using EarthSciDiscretizations

# Regenerate the rule catalog from `discretizations/` before each build so
# the manifest can never drift away from the JSON source of truth.
include(joinpath(@__DIR__, "generate_rule_catalog.jl"))
let repo_root = abspath(joinpath(@__DIR__, ".."))
    RuleCatalogGenerator.generate_rule_catalog(;
        catalog_dir = joinpath(repo_root, "discretizations"),
        output_path = joinpath(@__DIR__, "rule-catalog.md"),
        repo_root = repo_root,
    )
end

makedocs(
    sitename = "EarthSciDiscretizations.jl",
    modules = [EarthSciDiscretizations],
    authors = "EarthSciML Authors and Contributors",
    repo = "https://github.com/EarthSciML/EarthSciDiscretizations.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting started: solve a PDE" => "getting_started.md",
        "Finite-Volume Method" => "fv_method.md",
        "Cubed-Sphere Grid" => "grid.md",
        "Operators" => "operators.md",
        "Tutorial: Authoring a rule" => "tutorial.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://EarthSciML.github.io/EarthSciDiscretizations.jl",
        repolink = "https://github.com/EarthSciML/EarthSciDiscretizations.jl",
    ),
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/EarthSciML/EarthSciDiscretizations.jl",
    devbranch = "main",
)
