using AbstractMappings
using Documenter

DocMeta.setdocmeta!(AbstractMappings, :DocTestSetup, :(using AbstractMappings); recursive=true)

makedocs(;
    modules=[AbstractMappings],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="AbstractMappings.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/AbstractMappings.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/AbstractMappings.jl",
    devbranch="main",
)
