using Test, SafeTestsets

@safetestset "Linear algebra extension" begin

    using IgaBase, LinearAlgebra, SparseArrays

    @testset "Difference matrix" begin
        n = 4
        m = 3
        A = [-1  0  0;
              1 -1  0;
              0  1 -1;
              0  0  1]
        @test sparse_difference_matrix(n, Float64) == A
    end

    @testset "Nullspace of vector" begin

        @test permute_none([1,2,3,0,0,1,2,0,1,0]) == [1,2,3,6,7,9,4,5,8,10]

        # test zero vector
        a = zeros(10)
        C = nullspace(a)
        @test C==Matrix(I,10,10)

        # test with one non-zero entry
        a = zeros(10); a[4] = 2.0
        C = nullspace(a)
        @test size(C) == (10,9)
        @test isapprox(norm(a' * C), 0; atol=1e-12)

        # test dense vector
        a = rand(10) .- 0.5
        C = nullspace(a)
        @test size(C) == (10,9)
        @test isapprox(norm(a' * C), 0; atol=1e-12)

        # test sparse vector
        a = sprand(10, 0.5); a.nzval[:] .-= 0.5
        C = nullspace(a)
        @test size(C) == (10,9)
        @test isapprox(norm(a' * C), 0; atol=1e-12)
    end

    @testset "Nullspace of matrix" begin

        # test zero matrix
        A = zeros(10,10)
        C = nullspace(A)
        @test C==Matrix(I,10,10)

        # test dense matrix
        A = rand(5,10) .- 0.5
        C = nullspace(A)
        @test isapprox(norm(A * C), 0; atol=1e-12)

        # test sparse matrix
        A = sprand(5, 10, 0.5)
        C = nullspace(A)
        @test isapprox(norm(A * C), 0; atol=1e-12)
    end

end # safetestset
