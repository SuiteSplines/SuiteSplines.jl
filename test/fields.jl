# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module FieldsTest

using Test, LinearAlgebra
using AbstractMappings, StaticArrays, SortedSequences, CartesianProducts
using StaticArrays

@testset "Scalar Field" begin
    ϕ = Field((x,y,z) -> sin(π*x) * sin(π*y) * z)
    @test dimension(ϕ)==3
    @test codimension(ϕ)==(1,1)

    x = IncreasingRange(0.0,5.0,3) ⨱ IncreasingRange(0.0,2.0,4) ⨱ IncreasingRange(0.0,6.0,5)
    @evaluate y = ϕ(x)
    @test size(y)==size(x)
end

@testset "ScalarField iterable" begin
    ϕ = Field((x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)
    ϕ_test = [ϕ[1], ϕ[2], ϕ[3]]
    @test collect(ϕ) == ϕ_test
end

@testset "Pairing constructors" begin
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)
    alpha = Field{2,1}((x,y) -> x, (x,y) -> y)

    @test_throws AssertionError alpha ∘ F  
end

@testset "Trace of a vector valued field" begin
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)
    ϕ = Field((x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)

    X = IncreasingRange(0.0,2*π,4) ⨱ IncreasingRange(0.0,1.0,5)
    x = ∂(X,1)
    @evaluate z = (ϕ ∘ ∂(F, 1))(x)

    @test all(z.data[1] .≈ 2.0) && all(z.data[2] .≈ 0.0) && all(z.data[3] .≈ X.data[2])
end

@testset "Normal trace of a Gradient field" begin
    domain = Interval(0.0,1.0) ⨱ Interval(2.0,4.0)
    F = GeometricMapping(domain, (x,y) -> x, (x,y) -> y)
    ϕ = Field((x,y) -> x + y)

    # evaluation points on boundary 1
    X = CartesianProduct((d,n) -> IncreasingRange(d, n), domain, (3,4))
    x = ∂(X,1)

    # compute gradient restricted to boundary
    @evaluate y = (Gradient(ϕ) ∘ ∂(F, 1))(x)
    @test all(y.data[1] .≈ 1.0) && all(y.data[2] .≈ 1.0)

    # compute normal and normal trace of gradient of scalar field
    @evaluate n = Normal(∂(F, 1))(x)
    @test all(dot.(y, n) .≈ 1.0)
end

##
## POSSIBLY BROKEN?
##
#@testset "Pairing of a vector field and geometric mapping (Pullback)" begin
#    domain = Interval(0.0,2π) ⨱ Interval(0.0,1.0) 
#    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)
#    ϕ = Field{1, 2}((x,y,z) -> sin(π*x) * sin(π*y) * z, (x,y,z) -> sin(π*x) + cos(π*y)^2 * z)
#
#    x = CartesianProduct((d,n) -> IncreasingRange(d, n), domain, (3,4))
#    @evaluate y = F(x)
#    @evaluate z = ϕ(y)
#
#    @evaluate! z -= (ϕ ∘ F)(x)
#    @test norm(z) ≈ 0
#end

@testset "Gradient of a geometric mapping" begin
    domain = Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    F = Field{1,3}((x,y) -> x*y, (x,y) -> x*y + x + 1, (x,y) -> x*y + x + y + 2)
    G = Gradient(F)
    @test dimension(G) == 2
    @test codimension(G) == (2,3)

    x = IncreasingRange(0.0,2*π,3) ⨱ IncreasingRange(0.0,2.0,4)
    @evaluate y = G(x)
    @evaluate! y -= G(x)
    @test size(y) == size(x)
    @test sum(norm.(y)) ≈ 0
end

@testset "Hessian of a scalar mapping" begin
    domain = Interval(0.0,1.0) ⨱ Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    x = CartesianProduct((x,n) -> IncreasingRange(x,n), domain, (4,5,6))
    ϕ = Field((x,y,z) -> x^2 + y^2 + z^2)
    h = Hessian(ϕ)
    @evaluate y = h(x)
    @test y[2,2,2] == y[3,4,2] == Matrix(2*I,3,3)

    @evaluate! y -= h(x)
    @test size(y) == size(x)
    @test sum(norm.(y)) ≈ 0
end

end # FieldsTest