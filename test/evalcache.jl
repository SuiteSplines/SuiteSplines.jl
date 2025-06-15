using Test, SafeTestsets

@safetestset "Evaluation cache" begin

using TensorProductBsplines

space = SplineSpace(2,4) ⨷ SplineSpace(2,2) ⨷ SplineSpace(2,3);
spline = TensorProductBspline(space)
X = CartesianProduct(grevillepoints, spline.space)

import TensorProductBsplines: update!

@testset "update" begin
    Xᵣ = CartesianProduct(x -> global_insert(x,2), X)
    @test update!(spline.cache, Xᵣ) == true
    @test update!(spline.cache, Xᵣ) == false
    @test update!(spline.cache, X) == true
end

end
