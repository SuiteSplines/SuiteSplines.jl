export ⨱, CartesianProduct, collect!, dimension, dropdims

struct CartesianProduct{Dim,T,S<:Tuple} <: AbstractArray{T,Dim}
    data::S
    function CartesianProduct(v::AbstractVector)
        Dim = 1
        T = eltype(v)
        S = Tuple{typeof(v)}
        new{Dim,T,S}(tuple(v))
    end
    function CartesianProduct(args...)
        Dim = length(args)
        T = Tuple{map(eltype, args)...}
        S = typeof(args) 
        return new{Dim,T,S}(args)
    end
end
CartesianProduct(f::Type, args...) = CartesianProduct(map(f, args...)...)
CartesianProduct(f::Function, args...) = CartesianProduct(map(f, args...)...)

⨱(s...) = CartesianProduct(s...)
⨱(S₁::CartesianProduct, s₂) = CartesianProduct(S₁..., s₂)
⨱(s₁, S₂::CartesianProduct) = CartesianProduct(s₁, S₂...)
⨱(S₁::CartesianProduct, S₂::CartesianProduct) = CartesianProduct(S₁..., S₂...)

Base.convert(::Type{T}, x::T) where {T<:CartesianProduct{1}} = x
Base.convert(::Type{<:CartesianProduct{1}}, x) = CartesianProduct(x)

Base.ndims(::CartesianProduct{Dim}) where {Dim} = Dim 
Base.size(A::CartesianProduct) = ntuple(k -> length(A.data[k]), ndims(A))
Base.size(A::CartesianProduct, i) = length(A.data[i])
Base.eltype(::CartesianProduct{Dim,T}) where {Dim,T} = T

Base.iterate(X::CartesianProduct, state=1) = iterate(X.data, state)
Base.map(f, X::CartesianProduct) = Base.map(f, X.data)

@eval @inline function Base.getindex(X::CartesianProduct{1}, i::Int)
    return X.data[1][i]
end

for Dim in 2:4
    @eval @inline function Base.getindex(X::CartesianProduct{$Dim}, i::Vararg{Int, $Dim})
        return ntuple(k -> X.data[k][i[k]], $Dim)
    end
end

#Base.convert(::Type{S}, v::CartesianProduct{1,T,Tuple{S}}) where {T,S<:AbstractVector{T}} = v.data[1]  

using Base.Cartesian

"""
    collect(X)

Collect values of `CartesianProduct{T,Dim}` in `Dim`-dimensional
arrays. Similar to `ndgrid` in Matlab.

# Examples:
```jldoctest
julia> x, y = collect([1,2] ⨱ [2,3,4]);

julia> x
2×3 Matrix{Int64}:
 1  1  1
 2  2  2

julia> y
2×3 Matrix{Int64}:
 2  3  4
 2  3  4
```
"""
function Base.collect(X::CartesianProduct)
    Dim = ndims(X)
    A = map(x -> Array{eltype(x),Dim}(undef, size(X)), X)
    collect!(X, A)
    return A
end

@generated function collect!(X::CartesianProduct{Dim,T}, A::AbstractArray{T,Dim}) where {T,Dim}
    quote
        @assert length(X)==length(A)
        @nloops $Dim i A begin
            (@nref $Dim A i) = (@nref $Dim X i) 
        end
    end
end

@generated function collect!(X::CartesianProduct{Dim}, A::NTuple{Dim,<:AbstractArray}) where {Dim}
    quote
        @nexprs $Dim j->(A_j = A[j])
        @nexprs $Dim j->(size(A_j)==size(X))
        @nloops $Dim i A_1 d -> x_d = X.data[d][i_d] begin
            @nexprs $Dim j->(@nref $Dim A_j i) = x_j
        end
    end
end

# remove singular dimensions if there are any
function Base.dropdims(A::CartesianProduct; dims)
    @assert dims ∈ 1:ndims(A) && length(A.data[dims])==1
    return CartesianProduct(A.data[1:dims-1]..., A.data[dims+1:end]...)
end

function Base.convert(::Type{V}, X::CartesianProduct{1,T,Tuple{V}}) where {V,T}
    return X.data[1]
end

# get a specialized view on the data. Enables the @restrict macro.
@inline function IgaBase.restrict(X::CartesianProduct{Dim}, I::Vararg{Any,Dim}) where {Dim}
    return CartesianProduct((x,i) -> x[i], X.data, I)
end

