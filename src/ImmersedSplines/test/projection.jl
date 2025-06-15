module TestGalerkinProjection

using Test

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using TensorProductBsplines, IgaFormation, Algoim, ImmersedSplines
using LinearAlgebra, SparseArrays

@testset "Galerkin projection 2D" begin

    phi = AlgoimCallLevelSetFunction(
            (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2]]
        )

    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,2.0,11), IncreasingRange(0.0,2.0,11));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]);

    # project a polynomial onto the splinespace
    u_d = TensorProductBspline(U)
    f = ScalarFunction((x, y) -> 1 + x + y + x^2 + y^2 + y^3)
    project!(phi, u_d, onto=f, method=GalerkinProjection);

    # check polynomial reproduction
    @test l2_error(phi, u_d, to=f)[1] < 1e-12
end

@testset "Galerkin projection 3D" begin

    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2], 2.0 * x[3]]
        )

    # spline discretization
    partition = IncreasingRange(0.0,2.0,9) ⨱ IncreasingRange(0.0,2.0,10) ⨱ IncreasingRange(0.0,2.0,11);
    U = TensorProduct(u -> SplineSpace(3, u), partition);

    # project a polynomial onto the splinespace
    u_d = TensorProductBspline(U)
    f = ScalarFunction((x, y, z) -> 1 + x + y + 2z + x^2 + y^2 + 2z^2 + 0.5x^3 + y^3 + 0.5z^3)
    project!(phi, u_d, onto=f, method=GalerkinProjection);

    # check polynomial reproduction
    @test l2_error(phi, u_d, to=f)[1] < 1e-12
end


end