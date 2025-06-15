export order, get_matrices, get_row_sizes, get_col_sizes, AbstractKron
export issquare, issymmetric, isposdef

abstract type AbstractKron{T,N} <: AbstractArray{T,N} end

Base.eltype(::AbstractKron{T}) where T = T
Base.IndexStyle(::Type{<:AbstractKron}) = IndexCartesian()

"""
    order(M::AbstractMatrix)

Returns the order of a matrix, i.e. how many matrices are involved in the
Kronecker product (default to 1 for general matrices).
"""
IgaBase.order(::AbstractMatrix) = 1
IgaBase.order(A::AbstractKron) = length(A.data)

get_matrices(K::AbstractArray) = K
get_matrices(K::AbstractKron) = K.data

get_row_sizes(K::AbstractArray) = size(K,1)
get_row_sizes(K::AbstractKron)  = map(A -> size(A,1), get_matrices(K))
get_col_sizes(K::AbstractArray) = size(K,2)
get_col_sizes(K::AbstractKron)  = map(A -> size(A,2), get_matrices(K))
get_sizes(K) = get_row_sizes(K), get_col_sizes(K)

function get_kronecker_indices(i, remainder, nk)
    remainder = cld(remainder, nk)
    return cld(i,remainder), ((i-1) % remainder) + 1, remainder
end

"""
    Matrix(K::KroneckerProduct)

Converts a `GeneralizedKroneckerProduct` instance to a Matrix type.
"""
Base.Matrix(K::AbstractKron) = collect(K)

import Base: +

+(K1::AbstractKron, K2::AbstractKron) = collect(K1) + collect(K2)
function +(K::AbstractKron...)
    A = collect(K[1])
    for i in 2:length(K)
        A .+= collect(K[i])
    end
    return A
end

"""
    issquare(A::AbstractMatrix)

Checks if an array is a square matrix.
"""
function issquare(A::AbstractMatrix)
    m, n = size(A)
    return m == n
end

"""
    issymmetric(K::KroneckerProduct)

Checks if a Kronecker product is symmetric.
"""
function LinearAlgebra.issymmetric(K::AbstractKron)
    return check(issymmetric, K)
end

"""
    isposdef(K::KroneckerProduct)

Test whether a Kronecker product is positive definite (and Hermitian) by trying
to perform a Cholesky factorization of K.
"""
function LinearAlgebra.isposdef(K::AbstractKron)
    return check(isposdef, K)
end

function check(func, K::AbstractKron)
    return check(func, get_matrices(K))
end

# higher order function - used for simplifying checks on Kronecker matrices
function check(func, A)
    for a in A
        if !func(a)
            return false
        end
    end
    return true
end
