using Documenter
using TrixiAero

# Set up for local builds
if (get(ENV, "CI", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiAero, :DocTestSetup, :(using TrixiAero); recursive = true)

# Make documentation
makedocs(modules = [TrixiAero],
         sitename = "TrixiAero.jl",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                                  canonical = "https://acom-rwth.github.io/TrixiAero.jl/stable"),
         pages = [
             "Home" => "index.md",
             "API Reference" => "api.md"
         ],
         doctest = true,
         linkcheck = false,
         warnonly = [:missing_docs])

# Deploy documentation
deploydocs(repo = "github.com/ACoM-RWTH/TrixiAero.jl.git",
           devbranch = "main",
           push_preview = true)
