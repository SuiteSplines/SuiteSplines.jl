using Test

tests = [
    "base",
    "cartesian_product",
    "tensor_product"
]

@testset "CartesianProducts" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
