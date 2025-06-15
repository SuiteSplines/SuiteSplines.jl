using SortedSequences
using Documenter

DocMeta.setdocmeta!(SortedSequences, :DocTestSetup, :(using SortedSequences); recursive=true)

makedocs(;
    modules=[SortedSequences],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="SortedSequences.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/SortedSequences.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/SortedSequences.jl",
    devbranch="main",
)
