export Simplex, BezierSimplex
export b2c, c2b, barycentriccoords, cartesiancoords
export vol, jacobian, barycentricdir

"""
     Simplex{D,T<:Real} <: AbstractSimplex{D}

Simplex of dimension `D` with Vertices in real number type `T`.
"""
struct Simplex{D,T<:Real} <: AbstractSimplex{D}
    vertices::Frame{D,T}
end
Simplex(a::HPoint, b::HPoint...) = Simplex(Frame(a,b...))

_simplex(::Val{2}) = Simplex(HPoint(0.0), HPoint(1.0))
_simplex(::Val{3}) = Simplex(HPoint(0.0,0.0), HPoint(1.0,0.0), HPoint(0.0,1.0))
_simplex(::Val{4}) = Simplex(HPoint(0.0,0.0,0.0), HPoint(1.0,0.0,0.0), HPoint(0.0,1.0,0.0), HPoint(0.0,0.0,1.0))
Simplex(D::Dimension) = _simplex(Val(D))

"""
     BezierSimplex{D,T<:Real} <: AbstractSimplex{D}

Simplex of dimension `D` with Vertices in real number type `T`.
"""
struct BezierSimplex{P,D,T<:Real} <: AbstractSimplex{D}
    vertices::Frame{D,T}
    control::Vector{T}
end
BezierSimplex(P::Degree, f::Frame{D,T}, control::Vector{T}) where {D,T} = BezierSimplex{P,D,T}(f, control)
BezierSimplex(P::Degree, f::Frame{D,T}) where {D,T} = BezierSimplex(P, f, zeros(T, dimsplinespace(P,D)))
BezierSimplex(P::Degree, a::HPoint{D,T}, b::HPoint{D,T}...) where {D,T} = BezierSimplex(P, Frame(a,b...))

_beziersimplex(P, ::Val{2}) = BezierSimplex(P, HPoint(0.0), HPoint(1.0))
_beziersimplex(P, ::Val{3}) = BezierSimplex(P, HPoint(0.0,0.0), HPoint(1.0,0.0), HPoint(0.0,1.0))
_beziersimplex(P, ::Val{4}) = BezierSimplex(P, HPoint(0.0,0.0,0.0), HPoint(1.0,0.0,0.0), HPoint(0.0,1.0,0.0), HPoint(0.0,0.0,1.0))
BezierSimplex(P::Degree, D::Dimension) = _beziersimplex(P, Val(D))

"""
    dimension(simplex)

Return the dimension of a simplex.
"""
IgaBase.dimension(::AbstractSimplex{D}) where {D} = D

"""
    b2c(simplex)

Return the `mapping` (matrix) that maps barycentric coordinates to Cartesian coordinates.
"""
b2c(S::AbstractSimplex) = S.vertices
c2b(S::AbstractSimplex) = inv(b2c(S))

"""
    barycentriccoords(S::AbstractSimplex, x)

Compute the barycentric coordinates of a point, or matrix of points `x`.
"""
barycentriccoords(S::AbstractSimplex, x) = b2c(S) \ [x; ones(eltype(x),1,size(x,2))]

"""
    cartesiancoords(S::AbstractSimplex, λ)

Compute the Cartesian coordinates of a barycentric point or vector of barycentric points `λ`.
"""
cartesiancoords(S::AbstractSimplex, λ) = S.vertices * λ

"""
    vol(simplex)

Compute the volume of an affine simplex.
"""
vol(S::AbstractSimplex) = abs(det(b2c(S))) / factorial(dimension(S)-1)

"""
    jacobian(simplex)

Return the Jacobian mapping of a simplex.
"""
jacobian(S::AbstractSimplex) = c2b(S)[:,1:end-1]

barycentricdir(S::AbstractSimplex{D}, k=1:D-1) where {D} = c2b(S)[:,k]
