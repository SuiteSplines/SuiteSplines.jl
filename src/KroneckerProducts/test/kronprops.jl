@safetestset "KroneckerProduct base" begin

using KroneckerProducts, LinearAlgebra

@testset "properties" begin

    # associativity
    A, B, C = rand(2, 3), rand(3, 4), rand(4, 5)
    @test (A ⊗ B) ⊗ C == A ⊗ (B ⊗ C)

    # symmetry
    A = rand(4,4); A *=A'; B = rand(5,5); B *=B'
    @test issymmetric(A) && issymmetric(B)
    @test issymmetric(A ⊗ B)

    # positive definiteness
    @test isposdef(A) && isposdef(B)
    @test isposdef(A ⊗ B)

    # square
    @test issquare(A) && issquare(B)
    @test issquare(A ⊗ B)

    # other useful matrix properties
    @test det(A ⊗ B) ≈ det(Matrix(A ⊗ B))
    @test logdet(A ⊗ B) ≈ logdet(Matrix(A ⊗ B))
    @test tr(A ⊗ B) ≈ tr(Matrix(A ⊗ B))
end

@testset "Standard matrix operations" begin

    A, B = rand(2, 3), rand(3, 4)
    @test lmul!(2.5, copy(A) ⊗ B) ≈ 2.5 * Matrix(A ⊗ B)
    @test rmul!(A ⊗ copy(B), 2.5) ≈ 2.5 * Matrix(A ⊗ B)
    @test 2.5 * (A ⊗ B) ≈ 2.5 * Matrix(A ⊗ B)
    @test (A ⊗ B) * 2.5 ≈ 2.5 * Matrix(A ⊗ B)
end

@testset "Standard Matrix functions" begin

    A = rand(4,4); A *=A'; B = rand(5,5); B *=B'

    @test inv(A ⊗ B) == inv(A) ⊗ inv(B)
    @test transpose(A ⊗ B) == transpose(A) ⊗ transpose(B)
    @test adjoint(A ⊗ B) == adjoint(A) ⊗ adjoint(B)
    @test conj(A ⊗ B) == conj(A) ⊗ conj(B)

    A₁, B₁, C₁ = rand(2, 3), rand(3, 4), rand(4, 5)
    A₂, B₂, C₂ = rand(3, 2), rand(4, 3), rand(5, 4)

    # fast multiplication using mixed product property
    @test (A₁ ⊗ B₁ ⊗ C₁) * (A₂ ⊗ B₂ ⊗ C₂) ≈ Matrix(A₁ ⊗ B₁ ⊗ C₁) * Matrix(A₂ ⊗ B₂ ⊗ C₂)

    # fall-back method for standard matrix matrix multiplication 
    @test (A₁ ⊗ B₁) * (B₂ ⊗ A₂) ≈ Matrix(A₁ ⊗ B₁) * Matrix(B₂ ⊗ A₂)
end

end
