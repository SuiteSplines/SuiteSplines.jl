# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module BsplineRefinementTest

using Test

using LinearAlgebra
using TensorProductBsplines

@testset "TensorProductBspline patch refinement" begin
    space = SplineSpace(2, 4) ⨷ SplineSpace(2, 2) ⨷ SplineSpace(2, 3);
    spline = TensorProductBspline(space)
    spline.coeffs[:] = rand(size(spline.coeffs)...)

    spline₂ = refine(spline, method=kRefinement(2,2))
    x = CartesianProduct(breakpoints, spline.space)
    y = zeros(size(x))

    @evaluate! y = spline(x)
    @evaluate! y -= spline₂(x)

    @test isapprox(norm(y), 0.0, atol=1e-12)
end

end
