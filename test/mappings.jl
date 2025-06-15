# Normally we use SafeTestsets to create a safe scope for
# testing. However SafeTestsets cannot handle macros. So
# we use a module
module NurbsMappingTest

using Test
using SortedSequences, UnivariateSplines, TensorProductBsplines, NURBS, CartesianProducts
using LinearAlgebra: norm

@testset "Nurbs quarter circle" begin
    curve = arc(radius=1, origin=(0,0), α=0.0, β=2π/4)
    @test dimension(curve) == 1
    @test codimension(curve) == (1,2)

    x = IncreasingRange(0.0,1.0,5)
    @evaluate y = curve(x)
    radius² = norm.(y)
    @test all(radius² .≈ 1.0)
end

@testset "Three-dimensional Nurbs patch" begin
    space = SplineSpace(2,2) ⨷ SplineSpace(4, 3) ⨷ SplineSpace(3, 4);
    F = GeometricMapping(Nurbs, space; codimension=2)

    @test dimension(F) == 3
    @test codimension(F) == (1,2)

    # test partition of unity
    F.weights .= rand(size(F.space)...)
    F[1].coeffs .= F[2].coeffs .= 1.0
    transform_to_projective_space!(F)

    x = CartesianProduct(breakpoints, F.space)
    @evaluate y = F(x)
    @test norm(y.data[1] .- 1.0) ≈ 0 && norm(y.data[2] .- 1.0) ≈ 0
end


@testset "Evaluate on lines and surfaces of a 3D patch" begin

    mapping = cube(width=1, height=1, depth=1)
    Dim = 3
    surface = CartesianProduct([0.1,0.3],[0.1],[0.1,0.5])
    line = CartesianProduct([0.1],[0.1],[0.1,0.5])
    @evaluate f = mapping(line)
    @evaluate g = Gradient(mapping)(line)
    @evaluate f = mapping(surface)
    @evaluate g = Gradient(mapping)(surface)
end

@testset "Boundary Iterator Nurbs" begin
    space = SplineSpace(2,4) ⨷ SplineSpace(3,5) ⨷ SplineSpace(4,6)
    nurbs = Nurbs(space)

    # direct access to boundary components
    @test boundary(nurbs, 1) isa Nurbs{2}
    @test boundary(nurbs, 1,1) isa Nurbs{2}

    # boundary iterator
    B = Boundary(nurbs)
    @test eltype(B) <: Nurbs{2}
    @test length(B) == 6
    @test first(B) isa Nurbs{2}
end

end
