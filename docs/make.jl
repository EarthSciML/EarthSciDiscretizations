using Documenter
using EarthSciDiscretizations

makedocs(
    sitename = "EarthSciDiscretizations.jl",
    modules = Module[],
    pages = [
        "Home" => "index.md",
        "Finite-Volume Method" => "fv_method.md",
        "Cubed-Sphere Grid" => "grid.md",
        "Operators" => "operators.md",
        "Tutorial: Diffusion on the Sphere" => "tutorial.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
)
