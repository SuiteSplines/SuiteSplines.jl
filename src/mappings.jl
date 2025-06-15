export ScalarFunction, GeometricMapping, dimension, codimension

# function IgaBase.check_mapping_arguments(::Val{TensorProductBspline}, b::NTuple{Codim,TensorProductBspline}) where Codim
#     for k in 1:Codim
#         @assert sum(b[k].ders) == 0 "Provide B-spline functions; not their derivatives."
#     end
# end

# overloading of AbstractMappings.Evaluation.evalkernel! which enables computation
# on CartesianProduct grids of points using the @evaluation macro.
function AbstractMappings.evalkernel_imp!(op::Val{OP}, y, x, spline::AbstractSpline{Dim}) where {OP,Dim}
    update!(spline.cache, x)
    contract!(op, y, spline.coeffs, spline.cache.basis)
end

IgaBase.boundary_imp(spline::TensorProductBspline, comp, dir) = IgaBase.boundary_imp(spline, domain(spline), comp, dir)

# computation of boundary of a TensorProductBspline
function IgaBase.boundary_imp(spline::TensorProductBspline{2}, domain, comp, dir)
    space  = IgaBase.boundary_imp(spline.space,   comp, dir)
    ders   = IgaBase.boundary_imp(spline.ders,    comp, dir)
    coeffs = IgaBase.boundary_imp(spline.coeffs,  comp, dir)
    mindex = IgaBase.boundary_imp(spline.indices, comp, dir)
    return Bspline(space, coeffs, mindex; ders=ders[1], orientation=orientation(comp,dir))
end

function IgaBase.boundary_imp(spline::TensorProductBspline{3}, domain, comp, dir)
    space  = IgaBase.boundary_imp(spline.space,   comp, dir)
    ders   = IgaBase.boundary_imp(spline.ders,    comp, dir)
    coeffs = IgaBase.boundary_imp(spline.coeffs,  comp, dir)
    mindex = IgaBase.boundary_imp(spline.indices, comp, dir)
    return TensorProductBspline(space, coeffs, mindex; ders=ders, orientation=orientation(comp,dir))
end

function IgaBase.boundary_element_type(s::TensorProductBspline{2})
    J = Base.Slice{Base.OneTo{Int64}}
    T = numbertype(s.space)
    S = SubArray{T,1,typeof(s.coeffs),Tuple{J,Int64},true}
    I = SubArray{Int64,1,typeof(s.indices),Tuple{J,Int64},true}
    return Bspline{T,S,I}
end

function IgaBase.boundary_element_type(s::TensorProductBspline{3})
    J = Base.Slice{Base.OneTo{Int64}}
    T = numbertype(s.space)
    S = SubArray{T,2,typeof(s.coeffs),Tuple{J,J,Int64},true}
    I = SubArray{Int64,2,typeof(s.indices),Tuple{J,J,Int64},true}
    return TensorProductBspline{2,T,S,I}
end

IgaBase.num_boundary_components(::AbstractSpline{Dim}) where Dim = 2Dim
