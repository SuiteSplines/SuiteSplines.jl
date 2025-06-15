@safetestset "Box product identities" begin

using KroneckerProducts
using KroneckerProducts: box_exponent
using LinearAlgebra

@testset "Exponent in determinants" begin
    @test box_exponent(1,1) == 0
    @test box_exponent(2,1) == 0
    @test box_exponent(1,2) == 0
    @test box_exponent(2,2) == 1
end

@testset "Matrix operations" begin

    A, B = [1 2; 3 4], [5 6; 7 8]
    @test Matrix(A ⊠ B) == box(A,B)
    @test det(Matrix(A ⊠ B)) ≈ det(A ⊠ B)

    A = rand(2,2); A *=A'; B = rand(3,3); B *=B'
    @test inv(A ⊠ B) ≈ inv(Matrix(A ⊠ B))
    @test adjoint(A ⊠ B) ≈ adjoint(Matrix(A ⊠ B))
    @test transpose(A ⊠ B) ≈ transpose(Matrix(A ⊠ B))
    # @test conj(A ⊠ B) ≈ conj(Matrix(A ⊠ B)) (DOes not work yet)
end

A, B, C, D = rand(3,4), rand(4,5), rand(5,6), rand(4,5)

@testset "Box product identities" begin
    @test (A ⊠ B) ⊠ C == A ⊠ (B ⊠ C)
    @test (I(2) ⊠ I(3))' == I(3) ⊠ I(2)
    @test ((I(2) ⊠ I(3))') * (I(2) ⊠ I(3)) == I(6)
    @test (A ⊠ B) * (I(5) ⊠ I(4)) == A ⊗ B
    @test (A ⊗ B) * (I(4) ⊠ I(5)) == (A ⊠ B)
    @test (I(4) ⊠ I(3)) * (A ⊠ B) == B ⊗ A
    @test (I(4) ⊠ I(3)) * (A ⊗ B) == B ⊠ A
    @test (I(4) ⊠ I(3)) * vec(A') == vec(A)
    @test (I(3) ⊠ I(4)) * vec(A) == vec(A')
end

@testset "Mixed product properties" begin
    @test (A ⊠ B) * (C ⊠ D) ≈ Matrix(A ⊠ B) * Matrix(C ⊠ D)
    @test (A ⊗ B) * (D ⊠ C) ≈ Matrix(A ⊗ B) * Matrix(D ⊠ C)
    @test (A ⊠ B) * (C ⊗ D) ≈ Matrix(A ⊠ B) * Matrix(C ⊗ D)
end

end
