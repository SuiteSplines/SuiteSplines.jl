@safetestset "BoxProduct base" begin

using KroneckerProducts, LinearAlgebra

A, B, C = [1 2; 3 4], [5 6 7; 8 9 10], [11 12 13 14; 15 16 17 18; 19 20 21 22]

@testset "box-product" begin
    @test box(A) == A
    K = BoxProduct(A, A, B)
    @test size(K) == (8,12)
    @test order(K) == 3

    K = BoxProduct(A, B, C)
    @test order(K) == 3
    @test A ⊠ B ⊠ C == K
    @test A ⊠ B ⊠ C == collect(K) == Matrix(K) == box(A, B, C)

    @test BoxProduct(A,B,C; reverse=true) == C ⊠ B ⊠ A
    @test BoxProduct(n -> rand(n), 1:3) isa BoxProduct{Float64,1,3}
    @test BoxProduct(n -> rand(n,n+1), 1:3) isa BoxProduct{Float64,2,3}
    @test BoxProduct(n -> rand(n), 1:3; reverse=true) isa BoxProduct{Float64,1,3}
    @test BoxProduct(n -> rand(n,n+1), 1:3; reverse=true) isa BoxProduct{Float64,2,3}
    @test BoxProduct((n₁, n₂, n₃) -> rand(n₁, n₂) * rand(n₂, n₃), 1:3, 1:3, 1:3) isa BoxProduct{Float64,2,3}

end

@testset "getindex" begin
    K = box(A,B)
    C = [A[1,1]*B[1,1]  A[1,2]*B[1,1]   A[1,1]*B[1,2]   A[1,2]*B[1,2]     A[1,1]*B[1,3]   A[1,2]*B[1,3];
         A[1,1]*B[2,1]  A[1,2]*B[2,1]   A[1,1]*B[2,2]   A[1,2]*B[2,2]     A[1,1]*B[2,3]   A[1,2]*B[2,3];
         A[2,1]*B[1,1]  A[2,2]*B[1,1]   A[2,1]*B[1,2]   A[2,2]*B[1,2]     A[2,1]*B[1,3]   A[2,2]*B[1,3];
         A[2,1]*B[2,1]  A[2,2]*B[2,1]   A[2,1]*B[2,2]   A[2,2]*B[2,2]     A[2,1]*B[2,3]   A[2,2]*B[2,3]]
    @test K == C
end

end
