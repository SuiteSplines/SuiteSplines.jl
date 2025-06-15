# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module EvaluationTest

using Test, LinearAlgebra
using AbstractMappings, StaticArrays

@testset "base_2_ceil" begin
    import AbstractMappings.base_2_ceil
    @test base_2_ceil(Int, 127) == 128
    @test base_2_ceil(Int, 128) == 128
    @test base_2_ceil(Int, 129) == 256
    @test base_2_ceil(Int, 255) == 256
    @test base_2_ceil(Int, 256) == 256
    @test base_2_ceil(Int, 257) == 512
end

@testset "Grid evaluation 1d" begin

    import AbstractMappings.evalkernel!

    domain = Interval(0.0,3.0)
    X = IncreasingRange(domain, 4)
    f = x -> sin(x) * x^2 + x + 1
    f_at_x = f.(X)
    F = GeometricMapping(domain, f)

    @evaluate Y = F(X)
    @test Y==f_at_x

    @evaluate! Y += F(X)
    @test Y==2*f_at_x

    @evaluate! Y *= F(X)
    @test Y == 2 * (f_at_x.^2)

    @evaluate Y = F(X)
    @evaluate! Y -= F(X)
    @test Y == zeros(4)

    @evaluate Y = F(X)
    @evaluate! Y /= F(X)
    @test Y == ones(4)
end

@testset "Grid evaluation of a 3D geometry mapping" begin
    domain = Interval(0.0,3.0) ⨱ Interval(0.0,4.0) ⨱ Interval(0.0,5.0) 
    X = CartesianProduct((x,n) -> IncreasingRange(x, n), domain, (4,5,6))
    x,y,z = collect(X)
    h(x,y,z) = sin(x) * y^2 + z + 1
    h_at_x = h.(x,y,z)
    H = GeometricMapping(domain, h)

    @evaluate Y = H(X)
    @test Y == h_at_x

    @evaluate! Y += H(X)
    @test Y==2*h_at_x

    @evaluate! Y *= H(X)
    @test Y == 2*h_at_x.^2

    @evaluate Y = H(X)
    @evaluate! Y -= H(X)
    @test Y == zeros(size(X))

    @evaluate Y = H(X)
    @evaluate! Y /= H(X)
    @test Y == ones(size(X))
end

@testset "Grid evaluation of the gradient of a mapping" begin
    domain = Interval(0.0,3.0) ⨱ Interval(0.0,4.0) ⨱ Interval(0.0,3.0)
    x = CartesianProduct((x,n) -> IncreasingRange(x,n), domain, (4,5,6))
    F = GeometricMapping(domain, (x,y,z) -> sin(x) * y^2 + z + 1)
    G = Gradient(F)

    g = (x,y,z) -> [cos(x) * y^2, sin(x)*2*y,  1.0]

    @evaluate y = G(x)
    @test all(y[2,3,4] .≈ g(x[2,3,4]...))

    @evaluate! y -= G(x)
    @test norm(y.data[1]) + norm(y.data[2]) + norm(y.data[3]) ≈ 0
end

@testset "Grid evaluation of the hessian of a mapping" begin
    domain = Interval(0.0,3.0) ⨱ Interval(0.0,4.0) ⨱ Interval(0.0,5.0) 
    F = GeometricMapping(domain, (x,y,z) -> x^2 + y^2 + z^2 + 1)
    h = Hessian(F)

    x = CartesianProduct((x,n) -> IncreasingRange(x, n), domain, (4,5,6))
    @evaluate y = h(x)
    @test y[1] == y[10] == 2*Matrix(I,3,3)
end

end # EvaluationTest
