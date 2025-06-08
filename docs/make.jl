# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, GreenFunctionMonteCarlo

makedocs(
    modules = [GreenFunctionMonteCarlo],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Nils Niggemann",
    sitename = "GreenFunctionMonteCarlo.jl",

    pages = ["Overview" => "index.md","Contents" => "Contents.md","Tutorial" => "Example_transverseFieldIsing.md", "Reference.md"],
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl.git",
    push_preview = true
)
