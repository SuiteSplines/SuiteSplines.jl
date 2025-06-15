using Test

tests = [
    "core",
    "elementdata",
    "elementform",
    "patchdata",
    "patchform",
    "sumfact"
]

@testset "IgaFormation" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
