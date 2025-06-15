using Test

tests = [
    "algoim_interface",
    "elementform",
    "activefuns",
    "projection"
]

@testset "ImmersedSplines" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
