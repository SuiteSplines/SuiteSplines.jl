# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module NurbsProjectionTest

using Test
using SortedSequences, UnivariateSplines, TensorProductBsplines, NURBS, CartesianProducts
using LinearAlgebra: norm

dom = Interval(0.0,2.0) ⨱ Interval(0.0,3.0) ⨱ Interval(0.0,3.0)
space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,4,3), dom, (4,5,6));
nurbs = Nurbs(space)
nurbs.weights .= 2.0

f = ScalarFunction((x,y,z) -> 2*x^2 + 5*y^4*z^3)

@testset "Interpolation - polynomial reproduction" begin
    project!(nurbs, onto=f, method=Interpolation)
    x = CartesianProduct(breakpoints, nurbs.space)
    @evaluate y = nurbs(x)
    @evaluate! y -= f(x)
    @test isapprox(norm(y), 0.0, atol=1e-10)
end

@testset "Quasi interpolation - polynomial reproduction" begin
    project!(nurbs, onto=f, method=QuasiInterpolation)
    x = CartesianProduct(breakpoints, nurbs.space)
    @evaluate y = nurbs(x)
    @evaluate! y -= f(x)
    @test isapprox(norm(y), 0.0, atol=1e-10)
end

nurbs.weights .= rand(size(space)...) .+ 0.5

@testset "Interpolation - zero residual" begin
    project!(nurbs, onto=f, method=Interpolation)
    x = CartesianProduct(grevillepoints, nurbs.space)
    @evaluate y = nurbs(x)
    @evaluate! y -= f(x)
    @test isapprox(norm(y), 0.0, atol=1e-10)
end

f = GeometricMapping(dom, (x,y,z) -> 2*x^2 + 5*y^4*z^3, (x,y,z) -> 4*x^2*y^4*z^3)
nurbs = GeometricMapping(Nurbs, space; codimension=2)
nurbs.weights .= rand(size(space)...) .+ 0.5

@testset "Projection of mappings" begin
    project!(nurbs, onto=f, method=Interpolation)
    x = CartesianProduct(grevillepoints, nurbs.space)
    @evaluate y = nurbs(x)
    @evaluate! y -= f(x)
    @test isapprox(norm(norm.(y)), 0.0, atol=1e-10)
end

end
