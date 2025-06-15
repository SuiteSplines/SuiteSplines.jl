export refine, kRefinement, pRefinement, hRefinement, hpRefinement

function IgaBase.refine_imp(nurbs::Nurbs, method::AbstractRefinement)
    xw = IgaBase.refine_imp(nurbs.splinefun, method)
    w = nurbs.weightfun
    return Nurbs(xw, w)
end