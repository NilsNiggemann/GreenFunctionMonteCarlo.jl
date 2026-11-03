# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))

using Documenter, GreenFunctionMonteCarlo

makedocs(
    modules = [GreenFunctionMonteCarlo],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true",size_threshold_ignore = ["References.md"]),
    authors = "Nils Niggemann",
    sitename = "GreenFunctionMonteCarlo.jl",

    pages = ["Overview" => "index.md","Contents" => "Contents.md","Tutorial" => "Example_transverseFieldIsing.md", "Advanced Usage" => "Advanced.md","Reference" => "Reference.md"],
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
#
# Docs from `master` are deployed to `/dev`, docs from the `dev` branch to `/dev-branch`.
const IS_DEV_BRANCH = get(ENV, "GITHUB_REF", "") == "refs/heads/dev"

deploydocs(
    repo = "github.com/NilsNiggemann/GreenFunctionMonteCarlo.jl.git",
    devbranch = IS_DEV_BRANCH ? "dev" : "master",
    devurl = IS_DEV_BRANCH ? "dev-branch" : "dev",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev", "dev-branch" => "dev-branch"],
    push_preview = true
)
