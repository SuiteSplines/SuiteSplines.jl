export partial_derivative, Gradient, Hessian

@inline function IgaBase.partial_derivative(spline::TensorProductBspline{Dim}, ders::NTuple{Dim,Int}) where {Dim}
    return TensorProductBspline(spline.space, spline.coeffs, spline.indices; ders=spline.ders.+ders, orientation=spline.orientation)
end

@inline function IgaBase.partial_derivative(spline::TensorProductBspline{Dim}, k::Int, dir::Int) where Dim
    @assert (0 < dir < Dim+1)
    return partial_derivative(spline, ntuple(i -> k*Int(i==dir), Dim))
end