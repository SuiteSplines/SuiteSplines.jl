using Test

tests = [
            "base",
            "linalg",
            "quadrature"        
        ]

@testset "IgaBase" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
