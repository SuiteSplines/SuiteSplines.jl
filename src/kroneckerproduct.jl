export KroneckerProduct, ⊗, kron

"""
    KroneckerProduct{T,N,S}  <: AbstractVecOrMat{T}

Matrix type that stores the matrices in a Kronecker product and lazily
evaluates products, sums, etc.
"""
struct KroneckerProduct{T,N,M,S<:Tuple} <: AbstractKron{T,N}
    data::S
    sz::Tuple{Int,Int}
    function KroneckerProduct{N}(A::S) where {N, S<:Tuple}
        T = promote_type(map(eltype, A)...)
        M = length(A)
        m, n = 1, 1
        for a in A
            m *= size(a,1)
            n *= size(a,2)
        end
        return new{T,N,M,S}(A, (m, n))
    end
end

KroneckerProduct(A::AbstractVecOrMat...; reverse=false) = KroneckerProduct{2}((reverse==false) ? A : Base.reverse(A))
KroneckerProduct(A::AbstractVector...; reverse=false) = KroneckerProduct{1}((reverse==false) ? A : Base.reverse(A))
KroneckerProduct(f, iterable; reverse=false) = KroneckerProduct(map(f, iterable)...; reverse=reverse)
KroneckerProduct(f, iterable...; reverse=false) = KroneckerProduct(map(f, iterable...)...; reverse=reverse)

convert(::Type{KroneckerProduct{T,1,1}}, x::AbstractVector) where {T} = KroneckerProduct(x)

Base.size(K::KroneckerProduct) = K.sz
Base.size(K::KroneckerProduct, i) = K.sz[i]
Base.size(K::KroneckerProduct{T,1}) where {T} = (K.sz[1],)

⊗(A::AbstractVecOrMat, B::AbstractVecOrMat) = KroneckerProduct(A, B)
⊗(A::KroneckerProduct, B::AbstractVecOrMat) = KroneckerProduct(get_matrices(A)...,B)
⊗(A::AbstractVecOrMat, B::KroneckerProduct) = KroneckerProduct(A, get_matrices(B)...)

function Base.getindex(A::KroneckerProduct{T,1}, i1::Integer) where {T}
    m = size(A,1)
    value = T(1)
    for a in A.data
        ia, i1, m = get_kronecker_indices(i1, m, size(a,1))
        value *= a[ia]
    end
    return value
end

function Base.getindex(A::KroneckerProduct{T,2}, i1::Integer, i2::Integer) where {T}
    m, n = size(A)
    value = T(1)
    for a in A.data
        ia, i1, m = get_kronecker_indices(i1, m, size(a,1))
        ib, i2, n = get_kronecker_indices(i2, n, size(a,2))
        value *= a[ia, ib]
    end
    return value
end

using Base.Cartesian

# ToDo: there is quite some code repeated below. We can use code generation to
# generate the code for arbitrary dimension.

function Base.getindex(A::KroneckerProduct{T,2,2}, ::Colon, j::Integer) where T
    n = size(A,2)
    j_1, j, n = get_kronecker_indices(j, n, size(A.data[1],2))
    j_2, j, n = get_kronecker_indices(j, n, size(A.data[2],2))
    a_1 = getindex(A.data[1], :, j_1)
    a_2 = getindex(A.data[2], :, j_2)
    KroneckerProduct{1}((a_1, a_2)) # return view?
end

function Base.getindex(A::KroneckerProduct{T,2,3}, ::Colon, j::Integer) where T
    n = size(A,2)
    j_1, j, n = get_kronecker_indices(j, n, size(A.data[1],2))
    j_2, j, n = get_kronecker_indices(j, n, size(A.data[2],2))
    j_3, j, n = get_kronecker_indices(j, n, size(A.data[3],2))
    a_1 = getindex(A.data[1], :, j_1)
    a_2 = getindex(A.data[2], :, j_2)
    a_3 = getindex(A.data[3], :, j_3)
    KroneckerProduct{1}((a_1, a_2, a_3)) # return view?
end

function Base.getindex(A::KroneckerProduct{T,2,2}, i::Integer, ::Colon) where T
    m = size(A,1)
    i_1, i, m = get_kronecker_indices(i, m, size(A.data[1],1))
    i_2, i, m = get_kronecker_indices(i, m, size(A.data[2],1))
    a_1 = getindex(A.data[1], i_1, :)
    a_2 = getindex(A.data[2], i_2, :)
    KroneckerProduct{1}((a_1, a_2)) # return view?
end

function Base.getindex(A::KroneckerProduct{T,2,3}, i::Integer, ::Colon) where T
    m = size(A,1)
    i_1, i, m = get_kronecker_indices(i, m, size(A.data[1],1))
    i_2, i, m = get_kronecker_indices(i, m, size(A.data[2],1))
    i_3, i, m = get_kronecker_indices(i, m, size(A.data[3],1))
    a_1 = getindex(A.data[1], i_1, :)
    a_2 = getindex(A.data[2], i_2, :)
    a_3 = getindex(A.data[3], i_3, :)
    KroneckerProduct{1}((a_1, a_2, a_3)) # return view?
end

"""
    collect(K::KroneckerProduct)

Collects a lazy instance of the `KroneckerProduct` type into a full,
native matrix. Equivalent with `Matrix(K::KroneckerProduct)`.
"""
Base.kron(A::AbstractArray) = A
Base.collect(K::KroneckerProduct) = kron(get_matrices(K)...)