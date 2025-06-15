using Base.Cartesian

abstract type AbstractProduct{Dim} end

Base.ndims(A::AbstractProduct{Dim}) where {Dim} = Dim 
Base.size(A::AbstractProduct) = ntuple(k -> length(A[k]), ndims(A))
Base.size(A::AbstractProduct, i) = length(A.data[i])
Base.length(A::AbstractProduct) = ndims(A)

@inline function Base.getindex(X::AbstractProduct, i::Int)
    return X.data[i]
end

@inline Base.iterate(X::AbstractProduct, state=1) = iterate(X.data, state)

Base.map(f, X::AbstractProduct) = Base.map(f, X.data)

IgaBase.getdata(X::AbstractProduct, k) = X.data[k]