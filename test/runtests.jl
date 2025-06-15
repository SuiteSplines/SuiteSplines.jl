using Test
using SafeTestsets

tests = [
    "base",
    "kroneckerproduct",
    "kroneckersum",
    "kronprops",
    "factorizations",
    "contract",
    "boxproduct",
    "boxprops"
]

@testset "KroneckerProducts" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
