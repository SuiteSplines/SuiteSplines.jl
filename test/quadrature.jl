using Test, SafeTestsets

@safetestset "Quadrature" begin

    using IgaBase

    @testset "Construct quadrature rule" begin

        Q = QuadratureRule(rand(5), rand(5))
        @test Q isa AbstractQuadrule
        @test Q isa QuadratureRule
        @test length(Q) == 5
        @test size(Q) == (5,)
    end

end