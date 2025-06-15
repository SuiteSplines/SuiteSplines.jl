@safetestset "Base utilities" begin

using KroneckerProducts, LinearAlgebra

import KroneckerProducts: get_kronecker_indices

@testset "getindex" begin
    @test get_kronecker_indices(1, 20, 4) == (1, 1, 5)
    @test get_kronecker_indices(2, 20, 4) == (1, 2, 5)
    @test get_kronecker_indices(5, 20, 4) == (1, 5, 5)
    @test get_kronecker_indices(6, 20, 4) == (2, 1, 5)
    @test get_kronecker_indices(10, 20, 4) == (2, 5, 5)
    @test get_kronecker_indices(11, 20, 4) == (3, 1, 5)
end

@testset "Construction" begin

    A = rand(2,3) ⊗ rand(3,4)           # matrix
    B = rand(2) ⊗ rand(3) ⊗ rand(4)    # vector
    C = rand(3) ⊗ rand(4)'              # outer product

    @test A isa KroneckerProduct{Float64,2}
    @test B isa KroneckerProduct{Float64,1}
    @test C isa KroneckerProduct{Float64,2}

    @test order(A) == 2
    @test order(B) == 3
    @test order(C) == 2

    # inherited from AbstractArray
    @test ndims(A) == 2
    @test ndims(B) == 1
    @test ndims(C) == 2

    @test get_col_sizes(A) == (3,4)
    @test get_row_sizes(A) == (2,3)
    @test get_col_sizes(B) == (1,1,1)
    @test get_row_sizes(B) == (2,3,4)
    @test get_col_sizes(C) == (1,4)
    @test get_row_sizes(C) == (3,1)
end


@testset "addition" begin
    A = rand(2,3) ⊗ rand(3,4)
    B = rand(2,3) ⊗ rand(3,4)
    C = rand(2,3) ⊗ rand(3,4)
    @test A + B ≈ collect(A) + collect(B)
    @test A + B + C ≈ collect(A) + collect(B) + collect(C)
end

end
