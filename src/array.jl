export ResizableArray, maxlen, resize!

abstract type AbstractResizableArray{T,Dim} <:  DenseArray{T, Dim} end

Base.eltype(array::AbstractResizableArray) = eltype(array.data)
Base.ndims(array::AbstractResizableArray) = ndims(array.data)
Base.size(array::AbstractResizableArray) = size(array.data)
Base.length(array::AbstractResizableArray) = prod(size(array))

@inline Base.getindex(array::AbstractResizableArray, i...) = getindex(array.data, i...)
@inline Base.setindex!(array::AbstractResizableArray, v, i...) = setindex!(array.data, v, i...)

function Base.resize!(a::Array{T}, dims) where T
    @assert length(a) >= prod(dims)
    return unsafe_wrap(Array, pointer(a), dims)
end

"""
    ResizableArray{T,Dim} <: AbstractArray{T, Dim}

Array-type that may be resized. New resources will be allocated
if the new array-length exceeds maxlen. 
"""
mutable struct ResizableArray{T,Dim} <: AbstractResizableArray{T,Dim}
    raw::Vector{T}
    data::Array{T,Dim}
    function ResizableArray{T,Dim}(array_initializer, dims::NTuple{Dim,Int}) where {T,Dim}
        raw = Vector{T}(array_initializer, prod(dims))
        data = resize!(raw, dims)
        new{T,Dim}(raw, data)
    end
end

ResizableArray{T,Dim}(array_initializer, dims...) where {T,Dim} = ResizableArray{T,Dim}(array_initializer, dims) 
ResizableArray{T,Dim}(array_initializer; capacity::Int) where {T,Dim} = ResizableArray{T,Dim}(array_initializer, (capacity, ones(Int,Dim-1)...)) 

Base.pointer(array::ResizableArray) = pointer(array.data)
Base.pointer(array::ResizableArray, i::Integer) = pointer(array.data, i)
#Base.unsafe_convert(::Type{Ptr{T}}, array::ResizableArray{T}) where T = pointer(array.data)
#Base.elsize(array::ResizableArray) = elsize(array.data)

maxlen(array::ResizableArray) = length(array.raw)

function Base.resize!(array::ResizableArray{T, Dim}, dims::NTuple{Dim,Int}) where {T,Dim}
    @assert !(maxlen(array) < prod(dims)) # if the array is too small throw error
    array.data = unsafe_wrap(Array, pointer(array), dims)
end