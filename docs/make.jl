using Documenter
using EarthSciDiscretizations

makedocs(
    sitename = "EarthSciDiscretizations.jl",
    modules = [EarthSciDiscretizations],
    authors = "EarthSciML Authors and Contributors",
    repo = "https://github.com/EarthSciML/EarthSciDiscretizations.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
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
