export Interpolation, QuasiInterpolation, GalerkinProjection, project!, l2_error

import IgaBase.Interpolation, IgaBase.QuasiInterpolation, IgaBase.GalerkinProjection

function IgaBase.project!(spline::TensorProductBspline; onto, method::Type{<:AbstractInterpolation})
    x = CartesianProduct(grevillepoints, spline.space)
    update!(spline.cache, x)
    @evaluate! spline.coeffs = onto(x)
    IgaBase.project_imp!(method, spline, spline.coeffs)
    nothing
end

@inline function IgaBase.project_imp!(::Type{Interpolation}, spline::TensorProductBspline, y)
    @kronecker! y = inv(spline.cache.basis) * spline.coeffs
end

@inline function IgaBase.project_imp!(::Type{QuasiInterpolation}, spline::TensorProductBspline, y)
    A = KroneckerProduct(s -> Matrix(approximate_collocation_inverse(s)), spline.space; reverse=true)
    @kronecker! y = A * spline.coeffs
end