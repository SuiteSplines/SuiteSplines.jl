using NURBS
using Documenter

DocMeta.setdocmeta!(NURBS, :DocTestSetup, :(using NURBS); recursive=true)

makedocs(;
    modules=[NURBS],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="NURBS.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/NURBS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/NURBS.jl",
    devbranch="main",
)
