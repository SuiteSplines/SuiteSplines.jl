export TensorProductBspline, dimension, codimension, @evaluate!, orientation, grevillepoints
export domain, boundary, Boundary

IgaBase.numbertype(::TensorProduct{Dim,SplineSpace{T}}) where {Dim,T} = T

"""
    TensorProductBspline{Dim,Ders,T}(space, coeffs)
    TensorProductBspline{Dim,Ders,T}(space; codim=1)

Construct a `Dim` dimensional `TensorProductBspline` patch that is mapped to `Codim`
dimensional Euclidean space.
"""
struct TensorProductBspline{Dim,T<:Real,S<:AbstractArray{T,Dim}, I<:AbstractArray{Int,Dim}} <: AbstractSpline{Dim}
    space::TensorProduct{Dim,SplineSpace{T}}
    coeffs::S
    indices::I
    ders::NTuple{Dim,Int}
    cache::TPEvaluationCache{Dim}
    orientation::Int
    function TensorProductBspline(space::TensorProduct{Dim,SplineSpace{T}}, coeffs::S, indices::I, cache::TPEvaluationCache{Dim}; ders::NTuple{Dim,Int}, orientation::Int=1) where {Dim,T, S, I}
        @assert orientation==1 || orientation==-1
        @assert size(coeffs) == size(space) "Size coefficient array is inconsistent with dimension of spline-space."
        return new{Dim,T,S,I}(space, coeffs, indices, ders, cache, orientation)
    end
end

function TensorProductBspline(space::TensorProduct{Dim,SplineSpace{T}}, coeffs::AbstractArray{T,Dim}, indices::AbstractArray{Int, Dim}; ders=ntuple(k->zero(Int), Dim), orientation::Int=1) where {Dim,T}
    cache = TPEvaluationCache((s,d) -> (x -> Matrix(ders_bspline_interpolation_matrix(s, x, d+1)[d+1])), space, ders)
    return TensorProductBspline(space, coeffs, indices, cache; ders=ders, orientation=orientation)
end

function TensorProductBspline(space::TensorProduct{Dim,SplineSpace{T}}, coeffs::AbstractArray{T,Dim}, cache::TPEvaluationCache{Dim}; ders=ntuple(k->zero(Int), Dim), orientation::Int=1) where {Dim, T}
    indices = LinearIndices(map(s -> Base.OneTo(dimsplinespace(s)), space))
    return TensorProductBspline(space, coeffs, indices, cache; ders=ders, orientation=orientation)
end

function TensorProductBspline(space::TensorProduct{Dim,SplineSpace{T}}, coeffs::AbstractArray{T,Dim}; ders=ntuple(k -> zero(Int), Dim), orientation::Int=1) where {Dim, T}
    indices = LinearIndices(map(s -> Base.OneTo(dimsplinespace(s)), space))
    return TensorProductBspline(space, coeffs, indices; ders=ders, orientation=orientation)
end

function TensorProductBspline(space::TensorProduct{Dim,SplineSpace{T}}; orientation::Int=1) where {Dim, T}
    indices = LinearIndices(map(s -> Base.OneTo(dimsplinespace(s)), space))
    coeffs = zeros(T, size(space))
    return TensorProductBspline(space, coeffs, indices; orientation=orientation)
end

function Base.show(io::IO, spline::TensorProductBspline{Dim}) where Dim
    T = eltype(spline)
    S = spline.space
    println(io, "TensorProductBspline{$Dim} defined on splinespaces")
    for s in S
        println(io, "$s")
    end
end

function Base.similar(spline::TensorProductBspline; coeffs=similar(spline.coeffs))
    return TensorProductBspline(spline.space, coeffs, spline.indices, spline.cache; ders=spline.ders,orientation=spline.orientation)
end

IgaBase.orientation(spline::TensorProductBspline) = spline.orientation

IgaBase.domain(space::TensorProduct{Dim, <:SplineSpace}) where Dim = CartesianProduct(s -> domain(s), space)
IgaBase.domain(spline::TensorProductBspline) = domain(spline.space)

function UnivariateSplines.grevillepoints(space::TensorProduct{Dim,<:SplineSpace}) where Dim
    return CartesianProduct(grevillepoints, space)
end