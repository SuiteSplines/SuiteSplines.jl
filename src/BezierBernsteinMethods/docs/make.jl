using BezierBernsteinMethods
using Documenter

DocMeta.setdocmeta!(BezierBernsteinMethods, :DocTestSetup, :(using BezierBernsteinMethods); recursive=true)

makedocs(;
    modules=[BezierBernsteinMethods],
    authors="René Hiemstra, Michał Mika and contributors",
    sitename="BezierBernsteinMethods.jl",
    format=Documenter.HTML(;
        canonical="https://SuiteSplines.github.io/BezierBernsteinMethods.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SuiteSplines/BezierBernsteinMethods.jl",
    devbranch="main",
)
