using SuiteSplines
using Documenter

DocMeta.setdocmeta!(SuiteSplines, :DocTestSetup, :(using SuiteSplines); recursive=true)

makedocs(;
    modules=[SuiteSplines],
    authors="Michał Mika, René Hiemstra and contributors",
    sitename="SuiteSplines.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/SuiteSplines.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/SuiteSplines.jl",
    devbranch="main",
)
