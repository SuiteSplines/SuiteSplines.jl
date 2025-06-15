# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module BsplineTypesTest

using Test

using TensorProductBsplines

space = SplineSpace(2,4) ⨷ SplineSpace(2,2) ⨷ SplineSpace(3,3);
spline = TensorProductBspline(space)

@testset "Bspline patch contruction" begin
    @test dimension(spline) == 3
    @test codimension(spline) == (1,1)
end

x = CartesianProduct(grevillepoints, spline.space)
y = zeros(size(x))

@testset "Partition of unity" begin
    spline.coeffs .= 1.0
    @evaluate! y = spline(x)
    @test all(y.≈1.0)
end

@testset "Inplace operations" begin
    @evaluate! y += spline(x)
    @test all(isapprox.(y, 2.0, atol=1e-15))

    @evaluate! y -= spline(x)
    @test all(isapprox.(y, 1.0, atol=1e-15))
end

@testset "Boundary of a 2-dimension TensorProductBspline" begin
    s = TensorProductBspline(SplineSpace(2,5) ⨷ SplineSpace(3,6))
    s.coeffs[:] = 1:prod(size(s.coeffs))

    # direction 1, component 1
    ∂s = boundary(s, 1, 1)
    @test ∂s.coeffs == @view s.coeffs[1, :]
    @test ∂s.space == s.space[2]
    @test orientation(∂s) == -1

    # direction 1, component 2
    ∂s = boundary(s, 2, 1)
    @test ∂s.coeffs == @view s.coeffs[end, :]
    @test ∂s.space == s.space[2]
    @test orientation(∂s) == 1

    # direction 2, component 1
    ∂s = boundary(s, 1, 2)
    @test ∂s.coeffs == @view s.coeffs[:, 1]
    @test ∂s.space == s.space[1]
    @test orientation(∂s) == 1

    # direction 2, component 2
    ∂s = boundary(s, 2, 2)
    @test ∂s.coeffs == @view s.coeffs[:, end]
    @test ∂s.space == s.space[1]
    @test orientation(∂s) == -1
end

@testset "Boundary of a 3-dimensional TensorProductBspline" begin

    s = TensorProductBspline(SplineSpace(2,4) ⨷ SplineSpace(3,5) ⨷ SplineSpace(4,6))
    s.coeffs[:] = 1:prod(size(s.coeffs))

    # direction 1, component 1
    ∂s = boundary(s, 1, 1)
    @view s.coeffs[1, :, :]
    @test ∂s.coeffs == @view s.coeffs[1, :, :]
    @test ∂s.space[1] == s.space[2] && ∂s.space[2] == s.space[3]
    @test orientation(∂s) == -1

    # direction 1, component 2
    ∂s = boundary(s, 2, 1)
    @test ∂s.coeffs == @view s.coeffs[end, :, :]
    @test ∂s.space[1] == s.space[2] && ∂s.space[2] == s.space[3]
    @test orientation(∂s) == 1

    # direction 2, component 1
    ∂s = boundary(s, 1, 2)
    @test ∂s.coeffs == @view s.coeffs[:, 1, :]
    @test ∂s.space[1] == s.space[1] && ∂s.space[2] == s.space[3]
    @test orientation(∂s) == 1

    # direction 2, component 2
    ∂s = boundary(s, 2, 2)
    @test ∂s.coeffs == @view s.coeffs[:, end, :]
    @test ∂s.space[1] == s.space[1] && ∂s.space[2] == s.space[3]
    @test orientation(∂s) == -1

    # direction 2, component 1
    ∂s = boundary(s, 1, 3)
    @test ∂s.coeffs == @view s.coeffs[:, :, 1]
    @test ∂s.space[1] == s.space[1] && ∂s.space[2] == s.space[2]
    @test orientation(∂s) == -1

    # direction 2, component 2
    ∂s = boundary(s, 2, 3)
    @test ∂s.coeffs == @view s.coeffs[:, :, end]
    @test ∂s.space[1] == s.space[1] && ∂s.space[2] == s.space[2]
    @test orientation(∂s) == 1
end

import TensorProductBsplines: dir, comp

@testset "Boundary direction and component" begin
    @test dir(1) == 1 && comp(1) == 1
    @test dir(2) == 1 && comp(2) == 2
    @test dir(3) == 2 && comp(3) == 1
    @test dir(4) == 2 && comp(4) == 2
    @test dir(5) == 3 && comp(5) == 1
    @test dir(6) == 3 && comp(6) == 2
end

@testset "Boundary Iterator TensorProductBspline{2}" begin
    spline = TensorProductBspline(SplineSpace(3,5) ⨷ SplineSpace(4,6))
    a = boundary(spline, 4)
    b = boundary(spline, 2, 2)
    @test a.coeffs === b.coeffs && a.space === b.space

    B = Boundary(spline)
    @test eltype(B) == typeof(a)
    @test length(B) == 4
end

@testset "Boundary Iterator TensorProductBspline{3}" begin

    spline = TensorProductBspline(SplineSpace(2,4) ⨷ SplineSpace(3,5) ⨷ SplineSpace(4,6))
    a = boundary(spline, 5)
    b = boundary(spline, 1, 3)
    @test a.coeffs === b.coeffs && a.space === b.space

    B = Boundary(spline)
    @test eltype(B) == typeof(a)
    @test length(B) == 6
end

end
