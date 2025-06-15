using Test

tests = [
    "base",
    "mappings",
    "primitives",
    "refinement",
    "projection",
    "derivatives"
]

@testset "NURBS" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
