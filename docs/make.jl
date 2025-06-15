using CartesianProducts
using Documenter

DocMeta.setdocmeta!(CartesianProducts, :DocTestSetup, :(using CartesianProducts); recursive=true)

makedocs(;
    modules=[CartesianProducts],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="CartesianProducts.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/CartesianProducts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/CartesianProducts.jl",
    devbranch="main",
)
