using IgaBase
using Documenter

DocMeta.setdocmeta!(IgaBase, :DocTestSetup, :(using IgaBase); recursive=true)

makedocs(;
    modules=[IgaBase],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="IgaBase.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/IgaBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/IgaBase.jl",
    devbranch="main",
)
