@safetestset "KroneckerProduct base" begin

using KroneckerProducts, LinearAlgebra

@testset "Constructor" begin
    A, B = rand(4, 5), rand(5, 6)
    K = A ⊗ B
    @test size(K) == (size(A,1) * size(B,1), size(A,2) * size(B,2))
    @test order(K) == 2

    A, B, C = rand(2, 3), rand(3, 4), rand(4, 5)
    K = KroneckerProduct(A, B, C)
    @test order(K) == 3
    @test A ⊗ B ⊗ C == K
    @test A ⊗ B ⊗ C == collect(K) == Matrix(K) == kron(A, B, C)

    a, b, c = rand(2), rand(3), rand(4)
    K = KroneckerProduct(a, b, c)
    @test order(K) == 3
    @test a ⊗ b ⊗ c == K
    @test a ⊗ b ⊗ c == collect(K) == Matrix(K) == kron(a, b, c)
    @test KroneckerProduct(a,b,c; reverse=true) == c ⊗ b ⊗ a

    @test KroneckerProduct(n -> rand(n), 1:3) isa KroneckerProduct{Float64,1}
    @test KroneckerProduct(n -> rand(n,n+1), 1:3) isa KroneckerProduct{Float64,2}
    @test KroneckerProduct(n -> rand(n), 1:3; reverse=true) isa KroneckerProduct{Float64,1}
    @test KroneckerProduct(n -> rand(n,n+1), 1:3; reverse=true) isa KroneckerProduct{Float64,2}
    @test KroneckerProduct((n₁, n₂, n₃) -> rand(n₁, n₂) * rand(n₂, n₃), 1:3, 1:3, 1:3) isa KroneckerProduct{Float64,2}
end

@testset "getindex" begin
    A, B = rand(4, 5), rand(5, 6)
    K = KroneckerProduct(A, B);
    @test K[1,1] == A[1,1] * B[1,1]
    @test K[2,1] == A[1,1] * B[2,1]
    @test K[5,1] == A[1,1] * B[5,1]
    @test K[6,1] == A[2,1] * B[1,1]
    @test K[11,1] == A[3,1] * B[1,1]

    @test K[1,2] == A[1,1] * B[1,2]
    @test K[2,2] == A[1,1] * B[2,2]
    @test K[5,6] == A[1,1] * B[5,6]
    @test K[6,7] == A[2,2] * B[1,1]
    @test K[20,30] == A[4,5] * B[5,6]
    @test_throws BoundsError K[21,30]
    @test_throws BoundsError K[20,31]

    # get a column or row vector
    @test cat([K[:,k] for k in 1:size(K,2)]...; dims=2) == K
    @test cat([K[k,:] for k in 1:size(K,1)]...; dims=2)' == K

    A, B, C = rand(2, 3), rand(3, 4), rand(4, 5)
    K = KroneckerProduct(A, B, C)
    @test K[1,1] == A[1,1] * B[1,1] * C[1,1]
    @test K[4,1] == A[1,1] * B[1,1] * C[4,1]
    @test K[5,1] == A[1,1] * B[2,1] * C[1,1]
    @test K[8,1] == A[1,1] * B[2,1] * C[4,1]
    @test K[12,1] == A[1,1] * B[3,1] * C[4,1]
    @test K[13,1] == A[2,1] * B[1,1] * C[1,1]

    @test K[1,5] == A[1,1] * B[1,1] * C[1,5]
    @test K[1,6] == A[1,1] * B[1,2] * C[1,1]
    @test K[1,10] == A[1,1] * B[1,2] * C[1,5]
    @test K[1,15] == A[1,1] * B[1,3] * C[1,5]
    @test K[1,21] == A[1,2] * B[1,1] * C[1,1]

    @test_throws BoundsError K[25,1]
    @test_throws BoundsError K[1,61]

    # get a column or row vector
    @test cat([K[:,k] for k in 1:size(K,2)]...; dims=2) == K
    @test cat([K[k,:] for k in 1:size(K,1)]...; dims=2)' == K
end

end
