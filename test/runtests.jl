using Test

import SuiteSplines: SUITESPLINES_PKGS, bundle_test_include_mapexpr

@testset "SuiteSplines" begin
    for pkgname in SUITESPLINES_PKGS
        fp = joinpath(dirname(dirname(@__FILE__)), "src", pkgname, "test", "runtests.jl")
        println("$fp ...")
        include(bundle_test_include_mapexpr, fp)
    end
end # @testset
