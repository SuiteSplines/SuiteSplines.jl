using SuiteSplines
using Documenter

DocMeta.setdocmeta!(SuiteSplines, :DocTestSetup, :(using SuiteSplines); recursive=true)

makedocs(;
    modules=[SuiteSplines],
    doctest = false,
    authors="Michał Mika, René Hiemstra and contributors",
    sitename="SuiteSplines.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/SuiteSplines.jl",
        edit_link="main",
        assets=String[],
        size_threshold=1_000_000, # 1MB
    ),
    pages=[
        "About" => "index.md",
        "SpecialSpaces.jl" => "SpecialSpaces.md",
        "Index" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/SuiteSplines.jl",
    devbranch="main",
)
