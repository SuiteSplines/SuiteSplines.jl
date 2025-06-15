# Normally we use SafeTestsets to create a safe scope for
# testing. However SafeTestsets cannot handle macros. So
# we use a module
module NurbsPrimitivesTest

using Test, LinearAlgebra
using NURBS, CartesianProducts, UnivariateSplines, SortedSequences

@testset "Reference arc" begin
    pts = NURBS.cpts_reference_arc(1.0, π / 2)
    @test pts ≈ [1 1 0; 0 1 1]'

    wts = NURBS.wts_reference_arc(π / 2)
    @test wts ≈ [1, sqrt(2)/2, 1]

    @test NURBS.steps(π/4)==1
    @test NURBS.steps(2π/3)==1
    @test NURBS.steps(2.1π/3)==2
    @test NURBS.steps(4π/3)==2
    @test NURBS.steps(4.1π/3)==3
    @test NURBS.steps(2π)==3
end

@testset "Quarter arc" begin
    curve = arc(radius=1, origin=(0,0), α=0.0, β=2π/4)
    @test curve.weights ≈ [1, sqrt(2)/2, 1]
    x = grevillepoints(curve.space)

    @evaluate y = curve(x)
    @test all(norm.(y) .≈ 1)
end

@testset "One third arc" begin
    curve = arc(radius=1, origin=(0,0), α=0.0, β=2π/3)
    @test curve.weights ≈ [1, 1/2, 1]
    x = grevillepoints(curve.space)
    @evaluate y = curve(x)
    @test all(norm.(y) .≈ 1)
end

@testset "Circle" begin
    curve = circle(radius=1, origin=(0,0))
    x = grevillepoints(curve.space)
    @evaluate y = curve(x)
    @test all(norm.(y) .≈ 1)
end

curve = circle(radius=1, origin=(0,0))
s = curve.space

@testset "Annulus" begin
    o, rᵢ, rₒ = (1.0,2.0), 1.0, 2.0
    s = annulus(origin=o, inner_radius=rᵢ, outer_radius=rₒ)
    x = CartesianProduct(s -> grevillepoints(s), s[1].space)
    @evaluate y = s(x)

    z = sqrt.((y.data[1].-o[1]).^2 .+ (y.data[2].-o[2]).^2)
    @test isapprox(norm(z[:,1].-rₒ), 0.0, atol=1e-15)
    @test isapprox(norm(z[:,2].-rᵢ), 0.0, atol=1e-15)
end

@testset "Cylinder" begin
    o = (1.0,2.0,3.0)
    s = cylinder(origin=o, radius=1.0, height=2.0)
    x = CartesianProduct(grevillepoints, s.space)
    @evaluate y = s(x)
    @test all(sqrt.((y.data[1].-o[1]).^2 .+ (y.data[2].-o[2]).^2) .≈ 1.0)
end

@testset "Partial cylinder" begin
    o = (1.0,2.0,3.0)
    s = partial_cylinder(origin=o, radius=1.0, height=2.0, α=0.0, β=π/4)
    x = CartesianProduct(grevillepoints, s.space)
    @evaluate y = s(x)
    @test all(sqrt.((y.data[1].-o[1]).^2 .+ (y.data[2].-o[2]).^2) .≈ 1.0)
end

@testset "Tube" begin
    s = tube(inner_radius=0.25, outer_radius=0.75, height=1.0)
    x = CartesianProduct(grevillepoints, s.space)
    @evaluate y = s(x)
    @test all(sum(sqrt.(y.data[1].^2 .+ y.data[2].^2), dims=1) .≈ 1.0)
end

@testset "Partial tube" begin
    s = partial_tube(inner_radius=0.25, outer_radius=0.75, height=1.0, α = 0.0, β = π)
    x = CartesianProduct(grevillepoints, s.space)
    @evaluate y = s(x)
    @test all(sum(sqrt.(y.data[1].^2 .+ y.data[2].^2), dims=1) .≈ 1.0)
end

@testset "Rectangle" begin
    width, height = 10.0, 5.0
    s = rectangle(width=width, height=height)
    x = IncreasingVector([0, width]) ⨱ IncreasingVector([0, height])
    @evaluate y = s(x)

    @test s[1].coeffs == y.data[1]
    @test s[2].coeffs == y.data[2]
end

@testset "Cube" begin
    width, height, depth = 10.0, 5.0, 2.0
    v = cube(width=width, height=height, depth=depth)
    x = IncreasingVector([0, width]) ⨱ IncreasingVector([0, height]) ⨱ IncreasingVector([0, depth])
    @evaluate y = v(x)

    @test v[1].coeffs == y.data[1]
    @test v[2].coeffs == y.data[2]
    @test v[3].coeffs == y.data[3]
end

end
