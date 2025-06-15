# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module NurbsRefinementTest

using Test
using SortedSequences, UnivariateSplines, TensorProductBsplines, NURBS, CartesianProducts
using LinearAlgebra: norm

S = SplineSpace(2,2) ⨷ SplineSpace(4, 3) ⨷ SplineSpace(3, 4);
x = TensorProductBspline(S)
w = similar(x)
R = Nurbs(x,w)

@testset "Nurbs refinement" begin
    x.coeffs .=rand(size(S))
    w.coeffs .=rand(size(S))
    R₂ = refine(R, method=kRefinement(2,2))

    X = CartesianProduct(grevillepoints, R₂.space)
    @evaluate y = R(X)
    @evaluate y₂ = R₂(X)
    @test isapprox(norm(y-y₂), 0, atol=1e-12)
end

@testset "refinement of nurbs quarter circle" begin
    curve = arc(radius=1, origin=(0,0), α=0.0, β=2π/4)
    curve = refine(curve, method=kRefinement(2,2))
    x = grevillepoints(curve.space)
    @evaluate y = curve(x)
    radius² = norm.(y)
    @test all(radius² .≈ 1.0)
end

end
