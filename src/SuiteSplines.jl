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

for pkgname in SUITESPLINES_PKGS
    path = joinpath(pkgname, "src", "$pkgname.jl")
    include(bundle_include_mapexpr, path)
end

end # module SuiteSplines
