@safetestset "factorizations" begin

using KroneckerProducts, LinearAlgebra

A = rand(2,2); A *= A'
B = rand(3,3); B *= B'

@testset "LU" begin
    # LU currently does not work
end

@testset "Cholesky" begin
    K = A ⊗ B
    C = cholesky(K)
    C isa Cholesky{<:Real,<:KroneckerProduct}
    @test C.L * C.U ≈ K
    @test C.UL' * C.UL ≈ K
end

@testset "Eigen" begin
    K = A ⊗ B
    λ, ϕ = eigen(K)
    Λ, Φ = eigen(Matrix(K))

    @test ϕ * diagm(λ) * ϕ' ≈ Matrix(K)
    @test ϕ' * ϕ ≈ Matrix(I,6,6)
    @test sort(λ) ≈ sort(Λ)

    K = A ⊕ B
    λ, ϕ = eigen(K)
    Λ, Φ = eigen(Matrix(K))

    @test ϕ * diagm(λ) * ϕ' ≈ Matrix(K)
    @test ϕ' * ϕ ≈ Matrix(I,6,6)
    @test sort(λ) ≈ sort(Λ)
end

@testset "SVD" begin
    K = A ⊗ B
    U, S, V = svd(K)
    @test U * Diagonal(S) * V' ≈ K
    @test U' * U ≈ Matrix(I,6,6)
    @test V' * V ≈ Matrix(I,6,6)
end

end
