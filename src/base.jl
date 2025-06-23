export @suitesplines_reexport

"""
Vector of package names included in the `SuiteSplines.jl` bundle (without `.jl`).
"""
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
    "SpecialSpaces",
]

"""
    @suitesplines_reexport(pkgs...)

Convenience macro to reexport submodules of `SuiteSplines`.

- If no arguments are passed, it reexports all submodules in `SuiteSplines.SUITESPLINES_PKGS`.
- If one or more `Symbol`s are passed (e.g., `UnivariateSplines`), only those specific submodules are reexported.

This macro expands to a block of `@reexport using SuiteSplines.PkgName` statements.
"""
macro suitesplines_reexport(pkgs...)
    if length(pkgs) == 0
        return Expr(:block, [
            :( @reexport using SuiteSplines.$(Symbol(pkg)) )
            for pkg in SuiteSplines.SUITESPLINES_PKGS
        ]...)
    else
        return Expr(:block, [
            :( @reexport using SuiteSplines.$(pkg) )
            for pkg in pkgs
        ]...)
    end
end

function bundle_include_mapexpr(expr::Expr)
    MacroTools.postwalk(x ->
        begin
            !isa(x, Expr) && return x
            # if expression is include(path), make include(bundle_include_mapexpr, path)
            if x.head == :call
                if x.args[1] == :include
                    return :(include($bundle_include_mapexpr, $(x.args[2])))
                end
            end
            #@show x
            # clean up these loops...
            for pkgname in SUITESPLINES_PKGS
                # use modules in outer scope instead of dependencies
                x == :($(Expr(:., Symbol(pkgname)))) && return :($(Expr(:., :., :., Symbol(pkgname)))) 
                if x.head == :import
                    for k in eachindex(x.args)
                        if (x.args[k].args[1] == Symbol(pkgname))
                            x.args[k] = :($(Expr(:., :., :., Symbol(pkgname), x.args[k].args[2:end]...)))
                        end
                    end
                end
            end
            return x
        end,
    expr)
end

function bundle_test_include_mapexpr(expr::Expr)
    MacroTools.postwalk(x ->
        begin
            !isa(x, Expr) && return x
            # if expression is include(path), make include(bundle_include_mapexpr, path)
            if x.head == :call
                if x.args[1] == :include
                    return :(include($bundle_test_include_mapexpr, $(x.args[2])))
                end
            end
            # clean up these loops...
            for pkgname in SUITESPLINES_PKGS
                # use modules in outer scope instead of dependencies
                x == :($(Expr(:., Symbol(pkgname)))) && return :($(Expr(:., :SuiteSplines, Symbol(pkgname)))) 
                if x.head == :import
                    for k in eachindex(x.args)
                        if (x.args[k].args[1] == Symbol(pkgname))
                            x.args[k] = :($(Expr(:., :SuiteSplines, Symbol(pkgname), x.args[k].args[2:end]...)))
                        end
                    end
                end
            end
            return x
        end,
    expr)
end