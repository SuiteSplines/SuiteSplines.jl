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

function project_dependencies(project_toml_path::String)
    data = TOML.parsefile(project_toml_path)
    return get(data, "deps", Dict())
end

function collect_submodule_dependencies()
    deps = Dict{String,PackageSpec}()
    for pkg in SUITESPLINES_PKGS
        toml_path = joinpath("src", "bundle", pkg, "Project.toml")
        pkg_deps = project_dependencies(toml_path)
        for (name, uuid) in pkg_deps
            name in SUITESPLINES_PKGS && continue
            name in keys(deps) && continue
            deps[name] = PackageSpec(name=name, uuid=uuid) 
        end
    end
    return deps
end