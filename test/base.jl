# Normally we use SafeTestsets to create a safe scope for
# testing. However SafeTestsets cannot handle macros. So
# we use a module
module NurbsTypeTest

using Test
using LinearAlgebra
using NURBS, AbstractMappings, UnivariateSplines

import AbstractMappings.AbstractSpline

S = SplineSpace(2,2) ⨷ SplineSpace(4, 3) ⨷ SplineSpace(3, 4);
x = TensorProductBspline(S)
w = similar(x)
nurbs = Nurbs(x, w)
coeffs = ones(size(nurbs.space)...)
weights = rand(size(nurbs.space)...)

@testset "Nurbs construction" begin
    @test nurbs isa Nurbs{3}
    @test dimension(nurbs) == 3
    @test codimension(nurbs) == (1,1)
end

@testset "Nurbs partition of unity" begin
    nurbs.weights .= rand(size(nurbs.space)...)
    nurbs.coeffs .= 1.0
    transform_to_projective_space!(nurbs)
    x = CartesianProduct(breakpoints, nurbs.space)
    @evaluate y = nurbs(x)
    @test norm(y .- 1.0) ≈ 0
end

@testset "Projective space transformation" begin
    nurbs.weights .= weights
    nurbs.coeffs .= coeffs

    transform_to_projective_space!(nurbs)
    @test nurbs.coeffs == coeffs .* weights
    
    transform_from_projective_space!(nurbs)
    @test nurbs.coeffs == coeffs
end

@testset "Nurbs control points and weights" begin
    nurbs.weights .= weights
    nurbs.coeffs .= coeffs
    transform_to_projective_space!(nurbs)
    
    coeffs_test = UnivariateSplines.controlpoints(nurbs)
    weights_test = NURBS.nurbsweights(nurbs)

    @test coeffs_test ≈ coeffs
    @test weights_test ≈ weights
end

@testset "Nurbs refined control points and weights" begin
    nurbs.weights .= weights
    nurbs.coeffs .= coeffs
    transform_to_projective_space!(nurbs)

    R = refine(nurbs, method=kRefinement(2,2))
    coeffs_test = UnivariateSplines.controlpoints(R)
    @test all(coeffs_test .≈ 1.0)
end

@testset "Nurbs curve refined control points and weights" begin
    S = SplineSpace(2,2)
    x = Bspline(S)
    w = similar(x)
    nurbs = Nurbs(x, w)
    nurbs.coeffs .= ones(size(nurbs.space)...)
    nurbs.weights .= rand(size(nurbs.space)...)
    transform_to_projective_space!(nurbs)

    R = refine(nurbs, method=kRefinement(2,2))
    coeffs_test = UnivariateSplines.controlpoints(R)
    @test all(coeffs_test .≈ 1.0)
end

end
