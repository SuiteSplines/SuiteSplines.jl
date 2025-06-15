# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module MappingTest

using Test, LinearAlgebra
using AbstractMappings, StaticArrays, SortedSequences, CartesianProducts

@testset "Check number of input arguments" begin
    f(x,y,z) = x + 2y + 3z
    nargs = AbstractMappings.n_input_args(f)
    @test nargs isa Int
    @test nargs == 3

    g(x,y) = x + 2y
    nargs = AbstractMappings.n_input_args(g)
    @test nargs isa Int
    @test nargs == 2

    g(x) = x
    @test_throws AssertionError AbstractMappings.n_input_args(g)
end

@testset "Check domain type" begin
    @test AbstractMappings.isa_domain(Interval(0.0,1.0))
    @test AbstractMappings.isa_domain(Interval(0.0,1.0) ⨱ Interval(1.0,3.0))
end

@testset "GeometricMapping" begin
    # throw when dimensions are not consistent
    domain = Interval(0.0,1.0) ⨱ Interval(1.0,3.0)
    @test_throws ArgumentError GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (z) -> z)
end

@testset "GeometricMapping iterable" begin
    domain = Interval(0.0,1.0) ⨱ Interval(2.0,3.0) ⨱ Interval(4.0,5.0)
    F = GeometricMapping(domain, (x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)
    F_test = [F[1], F[2], F[3]]
    @test collect(F) == F_test
end

@testset "Cylinder" begin
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)
    @test dimension(F)==2
    @test codimension(F)==(1,3)

    x = IncreasingRange(0.0,2*π,10) ⨱ IncreasingRange(0.0,2.0,12) 
    @evaluate y = F(x)
    @test all(y.data[1].^2 + y.data[2].^2 .≈ 4)
end

@testset "boundary of a cylinder" begin
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)

    X = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (10,12))
    @evaluate Y = F(X)

    g = boundary(F, 1)
    x = X.data[2]
    @evaluate y = g(x)
    @test y[1] ≈ Y[1,1] && y[2] ≈ Y[1,2] && y[3] ≈ Y[1,3]

    g = boundary(F, 2)
    x = X.data[2]
    @evaluate y = g(x)
    @test y[1] ≈ Y[end,1] && y[2] ≈ Y[end,2] && y[3] ≈ Y[end,3]

    g = boundary(F, 3)
    x = X.data[1]
    @evaluate y = g(x)
    @test y[1] == Y[1,1] && y[2] == Y[2,1] && y[3] == Y[3,1]

    g = boundary(F, 4)
    x = X.data[1]
    @evaluate y = g(x)
    @test y[1] == Y[1,end] && y[2] == Y[2,end] && y[3] == Y[3,end]
end

@testset "boundary of a cube" begin
    domain = Interval(0.0,1.0) ⨱ Interval(2.0,3.0) ⨱ Interval(4.0,5.0)
    F = GeometricMapping(domain, (x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)

    X = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (10,11,12))
    @evaluate Y = F(X)

    g = boundary(F, 1)
    x = boundary(X,1)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[1,1,1] && y[end,end] ≈ Y[1,end,end]

    g = boundary(F, 2)
    x = boundary(X, 2)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[end,1,1] && y[end,end] ≈ Y[end,end,end]

    g = boundary(F, 3)
    x = boundary(X, 3)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[1,1,1] && y[end,end] ≈ Y[end,1,end]

    g = boundary(F, 4)
    x = boundary(X, 4)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[1,end,1] && y[end,end] ≈ Y[end,end,end]

    g = boundary(F, 5)
    x = boundary(X, 5)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[1,1,1] && y[end,end] ≈ Y[end,end,1]

    g = boundary(F, 6)
    x = boundary(X, 6)
    @evaluate y = g(x)
    @test y[1,1] ≈ Y[1,1,end] && y[end,end] ≈ Y[end,end,end]
end

@testset "Unit normal and volume element of curve" begin
    domain = Interval(0.0,2*π)
    F = GeometricMapping(domain, (θ) -> 2 * cos(θ), (θ) -> 2 * sin(θ))

    x = IncreasingRange(domain, 10) 
    @evaluate y = UnitNormal(F)(x)

    # test that normal is a unit-vector
    @test all(norm.(y) .≈ 1)

    # test for orthogonality w.r.t gradient vectors
    @evaluate ∇f = Gradient(F)(x)
    @test isapprox(norm(∇f .* y), 0.0, atol=1e-15)

    # evaluate the 'volume' element (Jacobian determinant) of the curve embedded in R²
    @evaluate y = Vol(F)(x)
    @test all(y .≈ 2)
end

@testset "Normal and volume element of cylinder" begin
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (θ,z) -> 2 * cos(θ), (θ,z) -> 2 * sin(θ), (θ,z) -> z)

    x = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (10,12)) 
    @evaluate y = UnitNormal(F)(x)

    # test that normal is a unit-vector
    @test all(norm.(y) .≈ 1)

    # test for orthogonality w.r.t gradient vectors
    @evaluate ∇f = Gradient(F)(x)
    @test isapprox(norm(∇f .* y), 0.0, atol=1e-15)

    # evaluate the 'volume' element (Jacobian determinant) of the surface embedded in R³    
    @evaluate y = Vol(F)(x)
    @test all(y .≈ 2)
end

@testset "Volume element of an annulus" begin
    domain = Interval(1.0,2.0) ⨱ Interval(0.0,2*π)
    F = GeometricMapping(domain, (r,θ) -> r * cos(θ), (r,θ) -> r * sin(θ))

    # evaluation points
    x = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (4,5)) 

    # evaluate the 'volume' element (Jacobian determinant) of the surface embedded in R³    
    @evaluate y = Vol(F)(x)
    @test all(y[1,:] .≈ 1)
    @test all(y[end,:] .≈ 2)
end

@testset "Volume element of a tube" begin
    domain = Interval(1.0,2.0) ⨱ Interval(0.0,2*π) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (r,θ,z) -> r * cos(θ), (r,θ,z) -> r * sin(θ), (r,θ,z) -> z)

    # evaluation points
    x = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (4,5,6)) 

    # evaluate the 'volume' element (Jacobian determinant) of the surface embedded in R³    
    @evaluate y = Vol(F)(x)
    @test all(y[1,:,:] .≈ 1)
    @test all(y[end,:,:] .≈ 2)
end

@testset "Gradient of a geometric mapping" begin
    domain = Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    F = GeometricMapping(domain, (x,y) -> x*y, (x,y) -> x*y + x + 1, (x,y) -> x*y + x + y + 2)
    G = Gradient(F)
    @test dimension(G) == 2
    @test codimension(G) == (2,3)

    x = IncreasingRange(0.0,2*π,10) ⨱ IncreasingRange(0.0,2.0,12)
    @evaluate y = G(x)
    @evaluate! y -= G(x)
    @test size(y) == size(x)
    @test sum(norm.(y)) ≈ 0
end

@testset "Hessian of a scalar mapping" begin
    domain = Interval(0.0,1.0) ⨱ Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    x = CartesianProduct((x,n) -> IncreasingRange(x,n), domain, (4,5,6))
    ϕ = GeometricMapping(domain, (x,y,z) -> x^2 + y^2 + z^2)
    h = Hessian(ϕ)
    @evaluate y = h(x)
    @test y[2,2,2] == y[3,4,2] == Matrix(2*I,3,3)

    @evaluate! y -= h(x)
    @test size(y) == size(x)
    @test sum(norm.(y)) ≈ 0
end

end # MappingTest