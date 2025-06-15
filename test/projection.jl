# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module BsplineProjectionTest

using Test

using LinearAlgebra
using TensorProductBsplines

space = SplineSpace(2, 4) ⨷ SplineSpace(1, 2) ⨷ SplineSpace(3, 3)
spline = TensorProductBspline(space)

@testset "Consistent interpolation" begin
    # interpolate a smooth function
    f = ScalarFunction((x,y,z) -> sin(x) * cos(2y) + z^2)
    project!(spline, onto=f, method=Interpolation)

    x = CartesianProduct(grevillepoints, spline.space)
    @evaluate y = f(x)
    @evaluate! y -= spline(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)
end

@testset "Consistent interpolation - polynomial reproduction" begin
    # reproduce a polynomial
    g = ScalarFunction((x,y,z) -> x^2 * y + z^3)
    project!(spline, onto=g, method=Interpolation)

    xᵣ = CartesianProduct(s -> global_insert(breakpoints(s), 2), spline.space)
    yᵣ = zeros(size(xᵣ))
    @evaluate! yᵣ = g(xᵣ)
    @evaluate! yᵣ -= spline(xᵣ)
    @test isapprox(norm(yᵣ), 0.0, atol=1e-12)

    @test isapprox(l2_error(spline, to=g)[1], 0.0, atol=1e-12)
    @test isapprox(l2_error(spline, to=g, relative=true)[1], 0.0, atol=1e-12)
end

@testset "Quasi interpolation - Lyche" begin
    # reproduce a polynomial
    g = ScalarFunction((x,y,z) -> x^2 * y + z^3)
    project!(spline, onto=g, method=QuasiInterpolation)

    x = CartesianProduct(s -> global_insert(breakpoints(s), 2), spline.space)
    y = zeros(size(x))
    @evaluate! y = g(x)
    @evaluate! y -= spline(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)

    # test L²error evaluation
    @test isapprox(l2_error(spline, to=g)[1], 0.0, atol=1e-12)
    @test isapprox(l2_error(spline, to=g, relative=true)[1], 0.0, atol=1e-12)
end

@testset "Project spline onto a spline" begin
    space₂ = TensorProduct(s -> SplineSpace(s.p, global_insert(s.U, 2)), spline.space);
    spline₂ = TensorProductBspline(space₂)    

    g = ScalarFunction((x,y,z) -> sin(x) * sin(y) * sin(z))
    project!(spline, onto=g, method=Interpolation)
    project!(spline₂, onto=spline, method=Interpolation)

    # check pointwise error
    x = CartesianProduct(s -> global_insert(breakpoints(s), 2), spline₂.space)
    y = zeros(size(x))
    @evaluate! y = spline(x)
    @evaluate! y -= spline₂(x)
    @test isapprox(norm(y), 0.0, atol=1e-12)

    # test L²error evaluation
    @test isapprox(l2_error(spline₂, to=spline)[1], 0.0, atol=1e-12)
    @test isapprox(l2_error(spline₂, to=spline, relative=true)[1], 0.0, atol=1e-12)
end

@testset "Project geometric mappings" begin
    dom = Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), dom, (4,5));

    g = GeometricMapping(dom, (u,v) -> 1 + u + u^2*v + u*v^3)
    gₕ = GeometricMapping(TensorProductBspline, space; codimension=1)

    # perform projection
    project!(gₕ, onto=g, method=Interpolation)

    # test L²error evaluation
    @test isapprox(l2_error(gₕ, to=g)[1], 0.0, atol=1e-12)
    @test isapprox(l2_error(gₕ, to=g, relative=true)[1], 0.0, atol=1e-12)
end

@testset "Project fields" begin
    dom = Interval(0.0,1.0) ⨱ Interval(0.0,1.0)
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), dom, (4,5));

    g = Field((u,v) -> 1 + u + u^2*v + u*v^3)
    gₕ = Field(TensorProductBspline, space)

    # perform projection
    project!(gₕ, onto=g, method=Interpolation)

    # test L²error evaluation
    @test isapprox(l2_error(gₕ, to=g)[1], 0.0, atol=1e-12)
    @test isapprox(l2_error(gₕ, to=g, relative=true)[1], 0.0, atol=1e-12)
end

end # module
