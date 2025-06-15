# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module BsplineMappingsTest

using Test

using TensorProductBsplines, LinearAlgebra

space = SplineSpace(2,4) ⨷ SplineSpace(2,2) ⨷ SplineSpace(3,3);
F = GeometricMapping(TensorProductBspline, space; codimension=2)

@testset "Mapping contruction" begin
    @test dimension(F) == 3
    @test codimension(F) == (1,2)
end

x = CartesianProduct(grevillepoints, space)

@testset "Inplace operations" begin
    F[1].coeffs .= F[2].coeffs .= 1.0
    @evaluate y = F(x)
    @test all(y[1].≈1.0) && all(y[2].≈1.0)

    @evaluate! y += F(x)
    @test all(isapprox.(y[1], 2.0, atol=1e-15))
    @test all(isapprox.(y[2], 2.0, atol=1e-15))

    @evaluate! y -= F(x)
    @test all(isapprox.(y[1], 1.0, atol=1e-15))
    @test all(isapprox.(y[2], 1.0, atol=1e-15))
end

@testset "Boundary 2D" begin

    # create analytical mapping
    domain = Interval(0.0,1.0) ⨱ Interval(2.0,3.0)
    f = GeometricMapping(domain, (x,y) -> x, (x,y) -> y)

    # create spline mapping
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), domain, (4,5));
    fₕ = GeometricMapping(TensorProductBspline, space; codimension=2)

    # project spline onto analytical mapping
    project!(fₕ, onto=f, method=Interpolation)

    # evaluation grid
    X = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (5,6))

    # check boundary implementation
    for k in 1:4
        k = 2
        x = boundary(X, k)
        b = boundary(f, k)
        bₕ = boundary(fₕ, k)

        @evaluate y = b(x)
        @evaluate! y -= bₕ(x)

        @test isapprox(norm(y), 0.0, atol=1e-12) 
    end
end

@testset "Boundary 3D" begin

    # create analytical mapping
    domain = Interval(0.0,1.0) ⨱ Interval(2.0,3.0) ⨱ Interval(4.0,5.0)
    f = GeometricMapping(domain, (x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)

    # create spline mapping
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3,4), domain, (4,5,6));
    fₕ = GeometricMapping(TensorProductBspline, space; codimension=3)

    # project spline onto analytical mapping
    project!(fₕ, onto=f, method=Interpolation)

    # evaluation grid
    X = CartesianProduct((d,n) -> IncreasingRange(d,n), domain, (10,11,12))

    # check boundary implementation
    for k in 1:6
        x = boundary(X, k)
        b = boundary(f, k)
        bₕ = boundary(fₕ, k)

        @evaluate y = b(x)
        @evaluate! y -= bₕ(x)

        @test isapprox(norm(y), 0.0, atol=1e-12) 
    end
end

end # module