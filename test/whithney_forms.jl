using Test, SafeTestsets

@safetestset "Whitney differential forms" begin

    using BezierBernsteinMethods
    import BezierBernsteinMethods: subsets, supp

    σ, τ = subsets((1,2,3), 1)
    @test length(σ) == length(τ) == binomial(3,1)
    @test σ[1] == (1,); @test τ[1] == (2,3)
    @test σ[2] == (2,); @test τ[2] == (1,3)
    @test σ[3] == (3,); @test τ[3] == (1,2)

    σ, τ = subsets((1,2,3,4), 2)
    @test length(σ) == length(τ) == binomial(4,2)
    @test σ[1] == (1,2); @test τ[1] == (3,4)
    @test σ[2] == (1,3); @test τ[2] == (2,4)
    @test σ[3] == (1,4); @test τ[3] == (2,3)
    @test σ[4] == (2,3); @test τ[4] == (1,4)
    @test σ[5] == (2,4); @test τ[5] == (1,3)
    @test σ[6] == (3,4); @test τ[6] == (1,2)

    # test condition in generating k-form patterns
    @test (0,1) ∪ supp([(1,0,0)...]) == [0,1]
    @test (0,1) ∪ supp([(0,1,0)...]) == [0,1]
    @test (0,1) ∪ supp([(0,0,1)...]) == [0,1,2]

    # test the pattern
    pattern = kformbasisconv(2, 1, 2)
    @test pattern.B == [1, 2, 2, 4, 3, 5, 1, 3, 2, 5, 3, 6, 4, 5, 5, 6]
    @test pattern.S ==  [1, -1/2, 1/2, -1, 1/2, -1/2, 1, -1/2, 1/2, -1/2, 1/2, -1, 1, -1/2, 1/2, -1]
    @test pattern.G[:] == [1,0,1,0,1,0,2,0,2,0,2,0,2,1,2,1]
end
