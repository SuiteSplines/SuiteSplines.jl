# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module BsplineDerivativesTest

using Test

using LinearAlgebra

using TensorProductBsplines

# construct spline
dom = Interval(0.0,4.0) ⨱ Interval(0.0,2.0) ⨱ Interval(0.0,3.0)
space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3,4), dom, (5,6,7));
spline = TensorProductBspline(space)

# coordinates to sample at
x = CartesianProduct(grevillepoints, space)
y = zeros(size(x))

@testset "derivatives" begin
    spline.coeffs .= 1.0

    # check derivatives of partition of unity are zero
    for k in 1:dimension(spline)
        ds = partial_derivative(spline, k, 1)
        @evaluate! y = ds(x)
        @test all(isapprox.(y, 0.0, atol=1e-15))
    end
end

# construct discrete and analytical mapping
gₕ = GeometricMapping(TensorProductBspline, space, codimension=1) # discrete solution
g = GeometricMapping(dom, (x,y,z) -> x^2 * y + z^3) # analytical solution

∇g = Gradient(g)   # analytical gradient
∇gₕ = Gradient(gₕ)  # discrete gradient

Δg = Hessian(g)    # analytical hessian
Δgₕ = Hessian(gₕ)   # discrete hessian

# perform projection and test accuracy w.r.t analytical solution
project!(gₕ, onto=g, method=Interpolation)

# construct B-spline scalar Field
@testset "gradient codim==1" begin

    # test gradient of Field
    @test ∇gₕ[1] isa TensorProductBspline{3} && ∇gₕ[1].ders == (1,0,0)
    @test ∇gₕ[2] isa TensorProductBspline{3} && ∇gₕ[2].ders == (0,1,0)
    @test ∇gₕ[3] isa TensorProductBspline{3} && ∇gₕ[3].ders == (0,0,1)

    # perform evaluation of complete Gradient
    @evaluate Z = ∇gₕ(x)
    @evaluate! Z -= ∇g(x)
    @test isapprox.(maximum(norm.(Z)), 0.0, atol=1e-12)
end

@testset "hessian" begin

    # test hessian of mapping
    @test Δgₕ[1,1] isa TensorProductBspline{3} && Δgₕ[1,1].ders == (2,0,0)
    @test Δgₕ[2,1] isa TensorProductBspline{3} && Δgₕ[2,1].ders == (1,1,0)
    @test Δgₕ[3,1] isa TensorProductBspline{3} && Δgₕ[3,1].ders == (1,0,1)

    @test Δgₕ[1,2] isa TensorProductBspline{3} && Δgₕ[1,2].ders == (1,1,0)
    @test Δgₕ[2,2] isa TensorProductBspline{3} && Δgₕ[2,2].ders == (0,2,0)
    @test Δgₕ[3,2] isa TensorProductBspline{3} && Δgₕ[3,2].ders == (0,1,1)

    @test Δgₕ[1,3] isa TensorProductBspline{3} && Δgₕ[1,3].ders == (1,0,1)
    @test Δgₕ[2,3] isa TensorProductBspline{3} && Δgₕ[2,3].ders == (0,1,1)
    @test Δgₕ[3,3] isa TensorProductBspline{3} && Δgₕ[3,3].ders == (0,0,2)

    # perform evaluation of complete Hessian
    @evaluate Z = Δgₕ(x)
    @evaluate! Z -= Δg(x)
    @test isapprox.(maximum(norm.(Z)), 0.0, atol=1e-10)
end

# @testset "gradient projection" begin
# # construct a 3D mapping
# ∇gₕ = GeometricMapping(TensorProductBspline, space; codimension=3)

# # project gradient of g
# project!(∇gₕ, onto=∇g, method=Interpolation)

# # compute jacobian and compare with hessian of g
# J = jacobian(∇s)

# # perform evaluation of complete Gradient
# Z = [zeros(size(X)) for i in 1:3, j in 1:3]
# @evaluate! Z = J(X)
# for j in 1:3
#     for i in 1:3
#         @evaluate! Y = Δg[i,j](X)
#         @test all(isapprox.(Y, Z[i,j], atol=1e-12))

#         @evaluate! Y -= J[i,j](X)
#         @test all(isapprox.(Y, 0.0, atol=1e-12))
#     end
# end
# end

end # end test-module
