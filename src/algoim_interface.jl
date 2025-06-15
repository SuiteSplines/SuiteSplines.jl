export ImmersedQuadRule

"""
    ImmersedQuadRule(map::AlgoimMapping, xa::Real, ya::Real, xb::Real, yb::Real, qo::Int64)

Compute a algoim quadrature rule in bounding box [xa, ya] × [xb, yb] based on a Gauss-Legendre
rule of `qo` points.
"""
struct ImmersedQuadRule{Dim, X, W}
    x::X
    w::W
    phase::Int
    function ImmersedQuadRule(phi::AlgoimCallLevelSetFunction, llcorner::NTuple{Dim,<:Real}, urcorner::NTuple{Dim,<:Real}; order::Int64, phase::Int64=-1) where {Dim}
        x, w = fill_quad_data(phi, SVector(llcorner), SVector(urcorner), phase, order)        
        X = typeof(x)
        W = typeof(w)
        return new{Dim,X,W}(x, w, phase) 
    end
end

ImmersedQuadRule(phi, element::Element{2}; order, phase=-1) = ImmersedQuadRule(phi, element[1,1], element[2,2]; order=order, phase=phase)
ImmersedQuadRule(phi, element::Element{3}; order, phase=-1) = ImmersedQuadRule(phi, element[1,1,1], element[2,2,2]; order=order, phase=phase)

# ToDo - this case is not implemented in AlgoimWrapper.jl
# function ImmersedQuadRule(phi, e::Element{1,<:SubArray}; order, phase=-1)
#     @show k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
#     s = e[1][k]
#     if k==1
#         f = (x) -> phi.φ(SVector(s,x[1]))
#         g = (x) -> phi.∇φ(SVector(s,x[1]))[1:1]
#         llcorner = (e[1][2],)
#         urcorner = (e[2][2],)
#     elseif k==2
#         f = (x) -> phi.φ(SVector(x[1],s))
#         g = (x) -> phi.∇φ(SVector(x[1],s))[1:1]
#         llcorner = (e[1][1],)
#         urcorner = (e[2][1],)
#     end
#     return ImmersedQuadRule(AlgoimCallLevelSetFunction(f,g), llcorner, urcorner; order=order, phase=phase)
# end

function ImmersedQuadRule(phi, e::Element{2,<:SubArray}; order, phase=-1)
    k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
    s = e[1,1][k]
    if k==1
        f = (x) -> phi.φ(SVector(s,x[1],x[2]))
        g = (x) -> phi.∇φ(SVector(s,x[1],x[2]))[1:2]
        llcorner = (e[1,1][2], e[1,1][3])
        urcorner = (e[2,2][2], e[2,2][3])
    elseif k==2
        f = (x) -> phi.φ(SVector(x[1],s,x[2]))
        g = (x) -> phi.∇φ(SVector(x[1],s,x[2]))[1:2]
        llcorner = (e[1,1][3], e[1,1][1])
        urcorner = (e[2,2][3], e[2,2][1])
    elseif k==3
        f = (x) -> phi.φ(SVector(x[1],x[2],s))
        g = (x) -> phi.∇φ(SVector(x[1],x[2],s))[1:2]
        llcorner = (e[1,1][1], e[1,1][2])
        urcorner = (e[2,2][1], e[2,2][2])
    end
    return ImmersedQuadRule(AlgoimCallLevelSetFunction(f,g), llcorner, urcorner; order=order, phase=phase)
end

Base.length(qr::ImmersedQuadRule) = Base.length(qr.w)
Base.size(qr::ImmersedQuadRule) = (length(qr),)