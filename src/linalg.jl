export sparse_difference_matrix, Permutation, PermuteNone, PermuteByAbsval
export nullspace, permute_by_abs, permute_none

using LinearAlgebra, SparseArrays

function sparse_difference_matrix(n::Int, T=Float64)
    return spdiagm(n, n-1, -1 => fill(T(1), n-1), 0 => fill(T(-1), n-1))
end

permute_by_abs(v::AbstractVector) = sortperm(v; by=abs)


function permute_none(v::AbstractVector{T}; atol=1e-13) where T
    m = length(v)
    indices = zeros(Int, m)
    a, b = 1, m
    for i in 1:m
        if abs(v[i])>atol
            indices[a] = i
            a+=1
        else
            indices[b] = i
            b-=1
        end
    end
    indices[b+1:end] = indices[end:-1:b+1]
    return indices
end


"""
    nullspace(A::AbstractMatrix{T<:Real})

Construct a left basis for the null space of matrix A.
"""
function LinearAlgebra.nullspace(A::AbstractMatrix{T}; perm=permute_by_abs) where {T<:Real}
    C = nullspace_matrix_imp(A; perm=perm)
    return C
end

function LinearAlgebra.nullspace(A::SparseMatrixCSC{T}; perm=permute_by_abs) where {T<:Real}
    C = nullspace_matrix_imp(A; perm=perm)
    droptol!(C, 1e2*eps(T))
    return C
end

function nullspace_matrix_imp(A::AbstractMatrix{T}; perm=permute_by_abs) where {T<:Real}
    m, n = size(A)

    # initialize induction step
    C = nullspace(A[1,:]; perm=perm)
    A = A * C

    for i in 2:m
        # compute basis for jth row of A
        temp = nullspace(A[i,:]; perm=perm)

        # update C
        C = C * temp

        # update A
        A = A * temp
    end
    return C
end

"""
    nullspace(a::AbstractVector{T}) where {T<:Real}

Find a sparse basis for the null space of a vector. The entries are computed
such that null space operator is maximally sparse and column sum of absolute
values equals 1.0 for improved conditioning.
"""
function LinearAlgebra.nullspace(a::SparseVector{T}; perm=permute_by_abs) where {T<:Real}
    return nullspace_impl(a; perm=perm)
end

function LinearAlgebra.nullspace(a::Vector{T}; perm=permute_by_abs) where {T<:Real}
    return Matrix(nullspace_impl(a; perm=perm))
end

# ToDo: Add option for two different update algorithms (defining c1)
function nullspace_impl(a::AbstractVector{T}; perm=permute_by_abs) where {T<:Real}

    # obtain permutation for more stable computations
    m = length(a)
    col = Vector{Int}(undef,2m)
    row = Vector{Int}(undef,2m)
    nzval = Vector{T}(undef,2m)
    p = perm(a)

    # treat initial zero values
    k, i, j  = 1, 1, 1
    while j < m+1 && isapprox(a[p[i]], 0; atol=100eps(T))
        row[k], col[k], nzval[k] = i, j, T(1)
        i+=1; j+=1; k+=1
    end

    # treat non-zero values
    c2 = T(0)
    while j < m && !isapprox(a[p[i+1]], 0; atol=100eps(T))
        f = a[p[i]] / a[p[i+1]]
        c1 = T(1) - c2
        # c1 = T(1) / (1. + abs(f))
        c2 = -c1 * f

        row[k],   col[k],   nzval[k]   = i,   j,  c1
        row[k+1], col[k+1], nzval[k+1] = i+1, j,  c2

        i += 1; j+=1; k += 2
    end

    # treat final zero values
    while i < m && isapprox(a[p[i+1]], 0; atol=100eps(T))
        row[k], col[k], nzval[k] = i+1, j, T(1)
        i += 1; j+=1; k += 1
    end

    # create and return sparse matrix
    n = (i>m) ? m : m-1
    return sparse(p[row[1:k-1]], col[1:k-1], nzval[1:k-1], m, n)
end
