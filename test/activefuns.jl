module TestExtendedSplines

using Test

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using TensorProductBsplines, IgaFormation, Algoim, ImmersedSplines
using LinearAlgebra, SparseArrays

@testset "Extension coefficiens p=2" begin
    # compute extension operator
    p = 2
    kts = KnotVector([0,0,0,1,3.0,3.0,3.0])
    n = dimsplinespace(p, kts)
    A = zeros(n, n-1)
    A[1,1] = A[2,2] = A[3,3] = 1
    A[4,:] = ImmersedSplines.compute_extension_coefficients(p, kts, 3, 4)

    # check for polynomial reproduction
    f = ScalarFunction((x) -> x^2 + x + 1)
    u = Bspline(SplineSpace(p, kts))
    project!(u, onto=f, method=Interpolation)
    y = (A' * A) \ (A' * u.coeffs)
    @test norm(A * y - u.coeffs[:]) < 1e-12
end

@testset "Extension coefficiens p=3" begin
    p = 3
    kts = KnotVector([0,0,0,0,1,3.0,3.0,3.0,3.0])
    n = dimsplinespace(p, kts)
    A = zeros(n, n-1)
    A[1,1] = A[2,2] = A[3,3] = 1
    A[4,:] = ImmersedSplines.compute_extension_coefficients(p, kts, 4, 4)
    A[5,:] = ImmersedSplines.compute_extension_coefficients(p, kts, 4, 5)

    # check for polynomial reproduction
    f = ScalarFunction((x) -> x^3 + x^2 + x + 1)
    u = Bspline(SplineSpace(p, kts))
    project!(u, onto=f, method=Interpolation)
    y = (A' * A) \ (A' * u.coeffs)
    @test norm(A * y - u.coeffs[:]) < 1e-12
end

@testset "Extension coefficiens p=2 with position constraint" begin
    # compute extension operator
    S = SplineSpace(2, IncreasingRange(0.0,3.0,4), cleft=1:1)
    n = dimsplinespace(S)
    A = zeros(n, n-1)
    A[2,1] = A[3,2] = A[4,3] = 1
    A[1,:] = ImmersedSplines.compute_extension_coefficients(S.p, S.U, 5, 2)
    A
    # check for polynomial reproduction
    f = ScalarFunction((x) -> x^2 + x)
    u = Bspline(S)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f, relative=true)[1] < 1e-12

    y = (A' * A) \ (A' * u.coeffs)
    @test norm(A * y - u.coeffs[:]) < 1e-12
end

@testset "Extension coefficiens p=2 with tangent constraint" begin
    # compute extension operator
    S = SplineSpace(2, IncreasingRange(0.0,3.0,4), cleft=2:2)
    n = dimsplinespace(S)
    A = zeros(n, n-1)
    A[2,1] = A[3,2] = A[4,3] = 1

    a = ImmersedSplines.compute_extension_coefficients(S.p, S.U, 5, 1)
    b = ImmersedSplines.compute_extension_coefficients(S.p, S.U, 5, 2)
    A[1,:] = (a+b)./2

    # check for polynomial reproduction
    f = ScalarFunction((x) -> x^2 + 1)
    u = Bspline(S)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f, relative=true)[1] < 1e-12

    y = (A' * A) \ (A' * u.coeffs)
    @test norm(A * y - u.coeffs[:]) < 1e-12
end

@testset "Extension operator 1D" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> -x[1]*x[1] + 1.5, 
        (x) -> -2.0 * x[1]
        )

    # spline discretization
    S = SplineSpace(2, IncreasingRange(0.0,3.0,4))
    C = spline_extension_operator(S, phi)

    # project a polynomial onto the splinespace
    u = Bspline(S)
    f = ScalarFunction((x) -> 1 + x + x^2)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f, relative=true)[1] < 1e-12

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs)) - u.coeffs
    indices = ImmersedSplines.active_splines(S, phi) .> 0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 1D with position constraint" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> -x[1]*x[1] + 1.5, 
        (x) -> -2.0 * x[1]
        )

    # spline discretization
    S = SplineSpace(2, IncreasingRange(0.0,3.0,4); cleft=1:1)
    C = spline_extension_operator(S, phi)

    # project a polynomial onto the splinespace
    u = Bspline(S)
    f = ScalarFunction((x) -> x + x^2)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f, relative=true)[1] < 1e-12

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs)) - u.coeffs
    indices = ImmersedSplines.active_splines(S, phi) .> 0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 1D with tangent constraint" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> -x[1]*x[1] + 1.5, 
        (x) -> -2.0 * x[1]
        )

    # spline discretization
    S = SplineSpace(2, IncreasingRange(0.0,3.0,4); cleft=2:2)
    C = spline_extension_operator(S, phi)

    # project a polynomial onto the splinespace
    u = Bspline(S)
    f = ScalarFunction((x) -> 1 + x^2)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f, relative=true)[1] < 1e-12

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs)) - u.coeffs
    indices = ImmersedSplines.active_splines(S, phi) .> 0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 2D" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,2.0,11), IncreasingRange(0.0,2.0,11));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]);
    C = spline_extension_operator(U, phi)
    
    # project a polynomial onto the splinespace
    u_d = TensorProductBspline(U)
    f = ScalarFunction((x, y) -> 1 + x + y + x^2 + y^2 + y^3)
    project!(u_d, onto=f, method=Interpolation)

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u_d.coeffs[:])) - u_d.coeffs[:]
    indices = ImmersedSplines.active_splines(U,phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 2D with position constraints" begin
    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.1,5), IncreasingRange(0.0,1.1,5));
    U = SplineSpace(2, partition.data[1]; cleft=1:1) ⨷ SplineSpace(3, partition.data[2]);
    C = spline_extension_operator(U, phi)

    
    # project a polynomial onto the splinespace
    u = TensorProductBspline(U)
    f = ScalarFunction((x, y) -> x + x*y + x^2 + x*y^2 + x*y^3)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f)[1] < 1e-14 # check that spline fits quadratic function exactly

    error = C * (C' * C \ (C' * u.coeffs[:])) - u.coeffs[:]
    indices = ImmersedSplines.active_splines(U,phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 2D with tangent constraints" begin
    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.1,5), IncreasingRange(0.0,1.1,5));
    U = SplineSpace(2, partition.data[1]; cleft=2:2) ⨷ SplineSpace(3, partition.data[2]; cleft=1:1);
    C = spline_extension_operator(U, phi)

    # project a polynomial onto the splinespace
    u = TensorProductBspline(U)
    f = ScalarFunction((x, y) -> (x^2+1) * (y + y^2 + y^3))
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f)[1] < 1e-14 # check that spline fits quadratic function exactly

    error = C * (C' * C \ (C' * u.coeffs[:])) - u.coeffs[:]
    indices = ImmersedSplines.active_splines(U,phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 3D" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2], 2.0 * x[3]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,2.0,11), IncreasingRange(0.0,2.0,11), IncreasingRange(0.0,2.0,11));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]) ⨷ SplineSpace(4, partition.data[3]);
    C = spline_extension_operator(U, phi)

    # project a polynomial onto the splinespace
    u = TensorProductBspline(U)
    f = ScalarFunction((x, y, z) -> 1 + x + y + z + x^2 + y^2 + z^2 + y^3 + z^3 + z^4)
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f)[1] < 1e-12 # check that spline fits quadratic function exactly

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs[:])) - u.coeffs[:]
    indices = ImmersedSplines.active_splines(U, phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 3D with mixed constraints" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2], 2.0 * x[3]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.1,7), IncreasingRange(0.0,1.1,8), IncreasingRange(0.0,1.1,9));
    U = SplineSpace(2, partition.data[1]; cleft=1:1) ⨷ SplineSpace(3, partition.data[2]; cleft=2:2) ⨷ SplineSpace(4, partition.data[3]; cleft=1:2);
    C = spline_extension_operator(U, phi)

    # project a polynomial onto the splinespace
    u = TensorProductBspline(U)
    f = ScalarFunction((x, y, z) -> (x + x^2) * (1 + y^2 + y^3) * (z^2 + z^3 + z^4))
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f)[1] < 1e-12 # check that spline fits quadratic function exactly

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs[:])) - u.coeffs[:]
    indices = ImmersedSplines.active_splines(U, phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end

@testset "Extension operator 3D with mixed constraints - spherical cavity" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2], 2.0 * x[3]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.1,7), IncreasingRange(0.0,1.1,8), IncreasingRange(0.0,1.1,9));
    U = TensorProduct(
        SplineSpace(2, partition.data[1]; cleft=1:1),
        SplineSpace(3, partition.data[2]),
        SplineSpace(4, partition.data[3])
    )
    C = spline_extension_operator(U, phi)

    # project a polynomial onto the splinespace
    u = TensorProductBspline(U)
    f = ScalarFunction((x, y, z) -> (x + x^2) * (1 + y + y^2 + y^3) * (1 + z + z^2 + z^3 + z^4))
    project!(u, onto=f, method=Interpolation)
    @test l2_error(u; to=f)[1] < 1e-12 # check that spline fits quadratic function exactly

    # check polynomial reproduction
    error = C * (C' * C \ (C' * u.coeffs[:])) - u.coeffs[:]
    indices = ImmersedSplines.active_splines(U, phi)[:].>0
    @test norm(error[indices])/sum(indices) < 1e-10
end


end