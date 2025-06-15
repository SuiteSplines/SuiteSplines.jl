using Test, SafeTestsets

@safetestset "Tensor product base" begin

using CartesianProducts

@testset "Tensor product constructors" begin
    s = (rand(3), rand(4), rand(5))
    S = TensorProduct(s...);

    @test ndims(S)==3
    @test S == TensorProduct(k -> s[k], 1:3)
    @test size(S) == (map(length, S)...,)
    @test S[1] == S.data[1] && S[2] == S.data[2] && S[3] == S.data[3]

    v = map(rand, 3:5)
    @test (v[1] ⨷ v[2]) ⨷ v[3] == v[1] ⨷ (v[2] ⨷ v[3])

    S = TensorProduct((a,b) -> a:b, 2:3, 4:5)
    @test S[1] == 2:4 && S[2] == 3:5
end

end
