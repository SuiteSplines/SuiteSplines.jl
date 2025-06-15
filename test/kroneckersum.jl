@safetestset "KroneckerSum" begin

using KroneckerProducts, LinearAlgebra

@testset "kronecker vectors sum" begin
    a, b = [1,2], [3,4,5]
    A, B = rand(3,3), rand(4,4)
    K = KroneckerSum(a,b)
    @test size(K) == (6,) && size(K,1) == 6 && size(K,2) == 1
    @test K == a ⊗ ones(Int,3) + ones(Int,2) ⊗ b
    @test K == a ⊕ b

    c = [6,7,8,9]
    K = KroneckerSum(a,b,c)
    @test size(K) == (24,) && size(K,1) == 24 && size(K,2) == 1
    @test K == a ⊗ ones(Int,3) ⊗ ones(Int,4) + ones(Int,2) ⊗ b  ⊗ ones(Int,4) + ones(Int,2) ⊗ ones(Int,3) ⊗ c
    @test K == a ⊕ b ⊕ c
    @test (a ⊕ b) ⊕ c == a ⊕ (b ⊕ c)
end

@testset "Kronecker matrix sum" begin

    A, B = [1 2; 3 4], [5 6 7; 8 9 10; 11 12 13]
    K = KroneckerSum(A, B)
    @test size(K) == (6,6) && size(K,1) == 6 && size(K,2) == 6
    @test get_row_sizes(K) == get_col_sizes(K) == (2,3)
    @test K == A ⊗ I(3) + I(2) ⊗ B
    @test A ⊕ B == K

    A, B, C = [1 2; 3 4], [5 6 7; 8 9 10; 11 12 13], [14 15 16 17; 18 19 20 21; 22 23 24 25; 26 27 28 29]
    K = KroneckerSum(A, B, C)
    @test size(K) == (24,24) && size(K,1) == 24 && size(K,2) == 24
    @test get_row_sizes(K) == get_col_sizes(K) == (2,3,4)
    @test K == A ⊗ I(3) ⊗ I(4) + I(2) ⊗ B ⊗ I(4) + I(2) ⊗ I(3) ⊗ C
    @test A ⊕ B ⊕ C == K
    @test (A ⊕ B) ⊕ C == A ⊕ (B ⊕ C)

end

end
