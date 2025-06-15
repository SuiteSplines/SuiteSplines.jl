using Test

tests = [
    "base",
    "evalcache",
    "types",
    "projection",
    "mappings",
    "derivatives",
    "refinement"
]

@testset "TensorProductBsplines" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
