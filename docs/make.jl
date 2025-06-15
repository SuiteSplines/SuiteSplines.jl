using IgaFormation
using Documenter

DocMeta.setdocmeta!(IgaFormation, :DocTestSetup, :(using IgaFormation); recursive=true)

makedocs(;
    modules=[IgaFormation],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="IgaFormation.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/IgaFormation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/IgaFormation.jl",
    devbranch="main",
)
