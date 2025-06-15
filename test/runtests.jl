using Test

tests = [
    "array_utils",
    "sequences"
]

if length(ARGS) > 0
    tests = ARGS
end

@testset "Sorted sequences" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
