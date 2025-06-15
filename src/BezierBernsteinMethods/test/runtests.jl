using Test

tests = [
    "bezier_types",
    "affine_geometry",
    "simplex",
    "multi_indexing",
    "bezierfuns",
    "bezier_methods",
    "whithney_forms"
]

@testset "BezierBernsteinMethods Unit-testing" begin

    # run all unit-tests
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end

end # @testset
