using Test

tests = [
    "array",
    "evalset",
    "evaluate",
    "mappings",
    "fields"
]

@testset "AbstractMappings" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
