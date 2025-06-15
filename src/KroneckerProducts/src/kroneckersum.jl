using LinearAlgebra: checksquare

export KroneckerSum, ⊕

struct KroneckerSum{T, N, M, S<:Tuple} <: AbstractKron{T,N}
    data::S
    sz::Int
    function KroneckerSum{1}(A::S) where {S<:Tuple} # useful for diagonal matrices
        T = promote_type(map(eltype, A)...)
        M = length(A)
        m = 1
        for a in A
            m *= length(a)
        end
        return new{T,1,M,S}(A, m)
    end
    function KroneckerSum{2}(A::S) where {S<:Tuple}
        T = promote_type(map(eltype, A)...)
        M = length(A)
        m = 1
        for a in A
            checksquare(a)
            m *= size(a,1)
        end
        return new{T,2,M,S}(A, m)
    end
end

@inline get_matrices(A::KroneckerSum) = A.data
@inline KroneckerSum(A::AbstractArray...) = KroneckerSum{2}(A)
@inline KroneckerSum(A::AbstractVector...) = KroneckerSum{1}(A)

⊕(A::AbstractArray, B::AbstractArray) = KroneckerSum(A, B)
⊕(A::KroneckerSum, B::AbstractArray)  = KroneckerSum(get_matrices(A)...,B)
⊕(A::AbstractArray, B::KroneckerSum)  = KroneckerSum(A, get_matrices(B)...)

Base.size(A::KroneckerSum{T,1}) where T = (A.sz,)
Base.size(A::KroneckerSum{T,2}) where T = (A.sz, A.sz)
Base.size(A::KroneckerSum{T,1}, i) where T = i>1 ? 1 : size(A)[i]
Base.size(A::KroneckerSum{T,2}, i) where T = size(A)[i]

@inline get_row_sizes(A::KroneckerSum) = map(a -> size(a,1), get_matrices(A))
@inline get_col_sizes(A::KroneckerSum) = map(a -> size(a,2), get_matrices(A))

function Base.getindex(A::KroneckerSum{T,1}, i1::Integer) where T
    m = length(A)
    value = T(0)
    for a in A.data
        ia, i1, m = get_kronecker_indices(i1, m, length(a))
        value += a[ia]
    end
    return value
end

function Base.getindex(A::KroneckerSum{T,2}, i1::Integer, i2::Integer) where T
    m, n = size(A)
    ia, ib, k = 1, 1, 1
    value = T(0)

    M = order(A)
    while ia==ib && k < M+1
        ia, i1, m = get_kronecker_indices(i1, m, size(A.data[k],1))
        ib, i2, n = get_kronecker_indices(i2, n, size(A.data[k],2))
        value += A.data[k][ia, ib] * (i1==i2)
        k+=1
    end
    return value
end
