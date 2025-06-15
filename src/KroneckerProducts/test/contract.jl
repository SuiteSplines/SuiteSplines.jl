# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module KroneckerContractionTests

using Test, KroneckerProducts

@testset "Matrix vector product in 1 dimension" begin
    K = KroneckerProduct(rand(5,5))
    X = rand(5)
    Y = zeros(5)
    Z = copy(X)

    contract!(Val(:(=)), Y, X, K)
    @test Y ≈ K * X

    # test in-place computation
    contract!(Val(:(=)), Z, Z, K)
    @test Z ≈ Y

    contract!(Val(:(+=)), Y, X, K)
    @test Y ≈ 2*Z

    contract!(Val(:(-=)), Y, X, K)
    @test Y ≈ Z
end

@testset "Kronecker dot product in 2 dimensions" begin
    n = (2,3)
    K = KroneckerProduct(map(rand, n)...)
    X = rand(n[end:-1:1]...)
    @test contract(X, K) ≈ (K.data[2]' * X * K.data[1])
    @test contract(X, K) ≈ dot(Vector(K), X[:])
end

@testset "Kronecker vector vector multiplication" begin
    m = (2,3,4,5)
    for k in 2:length(m)
        n = m[1:k]
        K = KroneckerProduct(map(rand, n)...)
        X = rand(n[end:-1:1]...)
        @test contract(X, K) ≈ dot(Vector(K), X[:])
        @test dot(K, X[:]) ≈ dot(Vector(K), X[:])
        @test dot(X[:], K) ≈ dot(Vector(K), X[:])
    end
end

@testset "Kronecker matrix vector multiplication 2d" begin
    E = zeros(2,3)
    X = rand(4,5)
    K = rand(3,5) ⊗ rand(2,4)

    @kronecker! E = K * X
    @test E ≈ K.data[2] * X * K.data[1]'
    @test E[:] ≈ Matrix(K) * X[:]

    @kronecker! E += K * X
    @test E ≈ 2 * K.data[2] * X * K.data[1]'
    @test E[:] ≈ 2 * Matrix(K) * X[:]

    @kronecker! E -= K * X
    @test E ≈ K.data[2] * X * K.data[1]'
    @test E[:] ≈ Matrix(K) * X[:]
end

end # module
