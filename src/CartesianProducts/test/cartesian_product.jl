using Test, SafeTestsets

@safetestset "Tensor product base" begin

using CartesianProducts

@testset "Constructor 1D" begin
    S = CartesianProduct(rand(4))
    @test S isa CartesianProduct{1,Float64,Tuple{Vector{Float64}}}
    @test ndims(S) == 1
    @test eltype(S) == Float64
    @test size(S) == (4,)
end

# #@testset "Conversion 1D" begin
# S = CartesianProduct(rand(4))
# convert(Vector{Float64}, S) == S.data[1]
# #end

@testset "Constructor" begin
    S = CartesianProduct(n -> rand(n), 2:4)
    @test S isa CartesianProduct{3,Tuple{Float64,Float64,Float64}}
    @test ndims(S) == 3
    @test eltype(S) == Tuple{Float64,Float64,Float64}
    @test size(S) == (2,3,4)
end

@testset "Cartesian base" begin
    u = map(k -> LinRange(0, k, k), 5:7)
    U = CartesianProduct(u[1], u[2], u[3])
    @test U == CartesianProduct(k -> LinRange(0, k, k), 5:7)

    @test size(U) == (5,6,7)
    @test length(U) == 210
    @test ndims(U) == 3
    @test eltype(U) == Tuple{Float64,Float64,Float64}
    @test (U.data[1] == u[1]) && (U.data[2] == u[2]) && (U.data[3] == u[3])
    @test U[end,end,end] == (5.0,6.0,7.0)

    @test u[1] ⨱ u[2] ⨱ u[3] == (u[1] ⨱ u[2]) ⨱ u[3] == u[1] ⨱ (u[2] ⨱ u[3])

    X = (1:4) ⨱ (2:5) ⨱ (4:9)
    @test X[1,2,3] == (1,3,6)

    for (k,x) in enumerate(X)
        @test x == X.data[k]
    end
end

@testset "Collect CartesianProduct 2D" begin
    X = CartesianProduct([1,2],[2,3,4])
    x, y = collect(X)  
    @test x == [1 1 1;
                2 2 2]
    @test y == [2 3 4;
                2 3 4]

    S = CartesianProduct((a,b) -> a:b, 2:3, 4:5)
    @test S[1,1] == (2,3)
    @test S[end,end] == (4,5)
end

@testset "Collect CartesianProduct 3D" begin
    X = CartesianProduct([1,2],[2,3,4], [5,6])
    @test ndims(X) == 3
    @test eltype(X) == Tuple{Int,Int,Int}
    @test size(X) == (2,3,2)
    @test size(X,1)==2 && size(X,2)==3 && size(X,3)==2
    @test X[end,end,end]==(2,4,6)

    S = CartesianProduct((a,b) -> a:b, 2:4, 4:6)
    @test S[1,1,1] == (2,3,4)
    @test S[end,end,end] == (4,5,6)
end

end
