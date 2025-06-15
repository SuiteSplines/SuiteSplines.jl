module SuiteSplines

using Pkg, TOML, Reexport, MacroTools

include("base.jl")

for pkgname in SUITESPLINES_PKGS
    path = joinpath(pkgname, "src", "$pkgname.jl")
    include(bundle_include_mapexpr, path)
end

end # module SuiteSplines
