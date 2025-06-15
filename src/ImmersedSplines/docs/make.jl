using ImmersedSplines
using Documenter

DocMeta.setdocmeta!(ImmersedSplines, :DocTestSetup, :(using ImmersedSplines); recursive=true)

makedocs(;
    modules=[ImmersedSplines],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="ImmersedSplines.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/ImmersedSplines.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/ImmersedSplines.jl",
    devbranch="main",
)
