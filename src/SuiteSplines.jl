module SuiteSplines

using Pkg, TOML, Reexport, MacroTools

const SUITESPLINES_PKGS = [
    "IgaBase",
    "SortedSequences",
    "CartesianProducts",
    "KroneckerProducts",
    "AbstractMappings",
    "BezierBernsteinMethods",
    "UnivariateSplines",
    "TensorProductBsplines",
    "NURBS",
    "IgaFormation",
    "ImmersedSplines",
]

include("base.jl")

if !haskey(ENV, "SUITESPLINES_PREPARE")
    for pkgname in SUITESPLINES_PKGS
        include(bundle_include_mapexpr, "$pkgname.jl")
    end
end

#@reexport using .IgaBase
#@reexport using .SortedSequences
#@reexport using .CartesianProducts
#@reexport using .KroneckerProducts
#@reexport using .AbstractMappings
#@reexport using .BezierBernsteinMethods
#@reexport using .UnivariateSplines
#@reexport using .TensorProductBsplines
#@reexport using .NURBS
#@reexport using .IgaFormation
#@reexport using .ImmersedSplines

end # module SuiteSplines
