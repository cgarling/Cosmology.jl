using Documenter
using Cosmology

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

DocMeta.setdocmeta!(Cosmology, :DocTestSetup, :(using Cosmology; import Unitful; import UnitfulAstro); recursive=true)

makedocs(
    sitename="Cosmology.jl",
    modules = [Cosmology],
    format = Documenter.HTML(;prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Chris Garling",
    pages = ["guide.md","types.md","public_methods.md","integrated_packages.md","private_methods.md","constants.md","index.md"],#,"api.md"],
    doctest=true
)

deploydocs(;
    repo = "github.com/cgarling/Cosmology.jl.git",
    versions = ["stable" => "v^", "v#.#"],
    push_preview=true,
)
