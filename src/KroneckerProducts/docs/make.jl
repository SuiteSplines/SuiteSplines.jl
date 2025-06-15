using KroneckerProducts
using Documenter

DocMeta.setdocmeta!(KroneckerProducts, :DocTestSetup, :(using KroneckerProducts); recursive=true)

makedocs(;
    modules=[KroneckerProducts],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="KroneckerProducts.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/KroneckerProducts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/KroneckerProducts.jl",
    devbranch="main",
)
