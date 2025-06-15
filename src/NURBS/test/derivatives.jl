# Normally we use SafeTestsets to create a safe scope for
# testing. However SafeTestsets cannot handle macros. So
# we use a module
module NurbsDerivativesTest

using Test
using LinearAlgebra: norm

using SortedSequences, UnivariateSplines, CartesianProducts, TensorProductBsplines
using NURBS

# construct a NURBS mapping
I = Interval(0.0,1.0)
S = SplineSpace(Degree(3), I, 3)
f = GeometricMapping(Nurbs, S ⨷ S ⨷ S; codimension=1)

# exact rational function of polynomials
s = ScalarFunction((x,y,z) ->  (x-1.0)^2*(y-0.2)*(z-2.0)^3)
w = ScalarFunction((x,y,z) ->  0.5 + 0.3x + 0.5y + 0.1z + 0.2x^2 + 0.4y^2 + 0.3z^2 + x*y*z)
g = GeometricMapping(I ⨱ I ⨱ I, (x,y,z) -> s(x,y,z) / w(x,y,z))

# Project Nurbs numerator and denominator
project!(f.splinefun, onto=s, method=Interpolation)
project!(f.weightfun, onto=w, method=Interpolation)

# test that rational is preserved
@test isapprox(l2_error(f, to=g)[1], 0.0, atol=1e-12)

# evaluation grid
x = CartesianProduct(breakpoints, f.space)

@testset "Nurbs evaluation" begin
    @evaluate y = f(x)
    @evaluate! y -= g(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)
end

@testset "Nurbs gradient evaluation" begin
    @evaluate y = Gradient(f)(x)
    @evaluate! y -= Gradient(g)(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)
end

@testset "Nurbs hessian evaluation" begin
    @evaluate y = Hessian(f)(x)
    @evaluate! y -= Hessian(g)(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)
end

@testset "Quarter circle tangent computation" begin
    f = arc(radius=1, origin=(0,0), α=0.0, β=2π/4)
    ∇f = Gradient(f)

    # compute position and tangent vector
    x = grevillepoints(f.space)
    @evaluate y = f(x)
    @evaluate ∇y = ∇f(x)

    # check orthogonality of radial vector and tangent vector
    @test isapprox(norm(y.data[1] .* ∇y.data[1] + y.data[2] .* ∇y.data[2]), 0.0, atol=1e-12)
end

end # end test-module
