using TensorProductBsplines
using Documenter

DocMeta.setdocmeta!(TensorProductBsplines, :DocTestSetup, :(using TensorProductBsplines); recursive=true)

makedocs(;
    modules=[TensorProductBsplines],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="TensorProductBsplines.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/TensorProductBsplines.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/TensorProductBsplines.jl",
    devbranch="main",
)
