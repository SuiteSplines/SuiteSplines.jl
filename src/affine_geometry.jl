export AffinePoint, AffineVector, AffinePointOrVector
export BPoint, BVector, HPoint, HVector, Frame

import Base: @propagate_inbounds, @boundscheck, checkbounds

abstract type AffinePoint{D,T} <: StaticVector{D, T} end
abstract type AffineVector{D,T} <: StaticVector{D, T} end
const AffinePointOrVector{D,T} = Union{AffinePoint{D,T}, AffineVector{D,T}}

struct BPoint{D,T<:Real} <: AffinePoint{D,T}
    data::NTuple{D,T}
    BPoint(X::NTuple{D,T}) where {D,T} = isapprox(sum(X), T(1); atol = eps(T)) ? new{D,T}(X) : throw(ArgumentError("Coordinates must sum to T(1)."))
end
BPoint(X::T...) where {T<:Real} = BPoint(X)

struct BVector{D,T<:Real} <: AffineVector{D,T}
    data::NTuple{D,T}
    BVector(X::NTuple{D,T}) where {D,T} = isapprox(sum(X), T(0); atol = eps(T)) ? new{D,T}(X) : throw(ArgumentError("Coordinates must sum to T(0)."))
end
BVector(X::T...) where {T<:Real} = BVector(X)

struct HPoint{D,T<:Real} <: AffinePoint{D,T}
    data::NTuple{D,T}
    HPoint(X::NTuple{D,T}) where {D,T} = isapprox(X[end], T(1); atol = eps(T)) ? new{D,T}(X) : throw(ArgumentError("Last coordinate must equal T(1)."))
end
HPoint(X::SVector) = HPoint(X.data)
HPoint(X::T...) where {T<:Real} = HPoint((X...,T(1)))

struct HVector{D,T<:Real} <: AffineVector{D,T}
    data::NTuple{D,T}
    HVector(X::NTuple{D,T}) where {D,T} = isapprox(X[end], T(0); atol = eps(T)) ? new{D,T}(X) : throw(ArgumentError("Last coordinate must equal T(0)."))
end
HVector(X::SVector) = HVector(X.data)
HVector(X::T...) where {T<:Real} = HVector((X...,T(0)))

@inline @propagate_inbounds function Base.getindex(A::AffinePointOrVector{D,T}, I::Vararg{Int,D}) where {D,T}
    @boundscheck checkbounds(A, I...)
    @inbounds ret = SVector(A.data[I...])
    return ret
end
@inline @propagate_inbounds function Base.getindex(A::AffinePointOrVector, I::AbstractUnitRange{Int})
    @boundscheck checkbounds(A, I)
    @inbounds ret = SVector(A.data[I])
    ret
end
@inline @propagate_inbounds function Base.getindex(A::AffinePointOrVector, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds ret = A.data[i]
    ret
end

import Base: -, +

(-)(p::BPoint{D,T}, q::BPoint{D,T}) where {D,T} = BVector(SVector(p)-SVector(q))
(-)(p::HPoint{D,T}, q::HPoint{D,T}) where {D,T} = HVector(SVector(p)-SVector(q))

(+)(p::BPoint{D,T}, v::BVector{D,T}) where {D,T} = BPoint(SVector(p)+SVector(v))
(+)(p::HPoint{D,T}, v::HVector{D,T}) where {D,T} = HPoint(SVector(p)+SVector(v))

const Frame{D,T<:Real} = SMatrix{D,D,T} # D-points forming a basis for D-dimensional space
Frame(a::HPoint{D,T}, b::HPoint{D,T}...) where {D,T} = Frame{D,T}([a b...])
