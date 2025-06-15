export Interpolation, QuasiInterpolation, project!

function IgaBase.project!(nurbs::Nurbs; onto, method::Type{<:AbstractInterpolation})
    x = grevillepoints(nurbs.space)
    transform_data_to_projective_space!(onto, nurbs, x)
    IgaBase.project_imp!(method, nurbs.splinefun, nurbs.coeffs)
end

function transform_data_to_projective_space!(f, nurbs::Nurbs, x)
    TensorProductBsplines.update!(nurbs.cache, x)
    @evaluate! nurbs.splinefun.coeffs = nurbs.weightfun(x)
    @evaluate! nurbs.splinefun.coeffs *= f(x)
end