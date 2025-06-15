export refine, hRefinement, pRefinement, hpRefinement, kRefinement

function IgaBase.refine_imp(spline::TensorProductBspline, method::AbstractRefinement)
    # determine refinement operator and new spline spaces
    extraction_operators, univariate_spaces = unzip(map(s -> refinement_operator(s, method), spline.space))
    C = KroneckerProduct(c -> Matrix(c), extraction_operators; reverse=true)
    space = TensorProduct(univariate_spaces...)

    # compute new coefficients
    coeffs = zeros(eltype(spline.coeffs), size(space))
    @kronecker! coeffs = C * spline.coeffs

    return TensorProductBspline(space, coeffs; ders=spline.ders, orientation=spline.orientation)
end
