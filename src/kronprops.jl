using LinearAlgebra: checksquare

export det, logdet, tr, lmul!, rmul!

"""
    det(K::KroneckerProduct)

Compute the determinant of a Kronecker product.
"""
function LinearAlgebra.det(K::KroneckerProduct{T,2}) where T
    checksquare(K)
    if check(issquare, get_matrices(K))
        m = size(K,1)
        determinant = one(T)
        for a in get_matrices(K)
            x = cld(m, size(a,1))
            determinant *= det(a)^x
        end
        return determinant
    else
        return zero(T)
    end
end

"""
    logdet(K::KroneckerProduct)

Compute the logarithm of the determinant of a Kronecker product.
"""
function LinearAlgebra.logdet(K::KroneckerProduct{T,2}) where T
    checksquare(K)
    if check(issquare, get_matrices(K))
        m = size(K,1)
        determinant = zero(T)
        for a in get_matrices(K)
            x = cld(m, size(a,1))
            determinant += x * logdet(a)
        end
        return determinant
    else
        return real(T)(-Inf)
    end
end

"""
    tr(K::KroneckerProduct)

Compute the trace of a Kronecker product.
"""
function LinearAlgebra.tr(K::KroneckerProduct{T,2}) where T
    checksquare(K)
    if check(issquare, get_matrices(K))
        trace = one(T)
        for A in get_matrices(K)
            trace *= tr(A)
        end
        return trace
    else
        return sum(diag(K))
    end
end

"""
    inv(K::KroneckerProduct)

Compute the inverse of a Kronecker product.
"""
function Base.inv(K::KroneckerProduct)
    checksquare(K)
    if check(issquare, get_matrices(K))
        return KroneckerProduct(map(inv, get_matrices(K))...)
    else
        throw(SingularException(1))
    end
end

"""
    adjoint(K::KroneckerProduct)

Compute the adjoint of a Kronecker product.
"""
function Base.adjoint(K::KroneckerProduct)
    return KroneckerProduct(map(adjoint, get_matrices(K))...)
end

"""
    transpose(K::KroneckerProduct)

Compute the transpose of a Kronecker product.
"""
function Base.transpose(K::KroneckerProduct)
    return KroneckerProduct(map(transpose, get_matrices(K))...)
end

function Base.conj(K::AbstractKron)
    return KroneckerProduct(map(conj, get_matrices(K))...)
end

function Base.:*(A::KroneckerProduct, B::KroneckerProduct)
    if get_col_sizes(A) == get_row_sizes(B)
        # mixed product poperty
        return KroneckerProduct{2}(ntuple(k -> A.data[k] * B.data[k], order(A)))
    end
    @warn "Mismatch between A and C in (A ⊗ B)(C ⊗ D), using less efficient fall-back routines."
    return Matrix(A) *  Matrix(B)
end

"""
    lmul!(a::Number, K::KroneckerProduct)

Scale a `KroneckerProduct` `K` inplace by a factor `a` by rescaling the
**left** matrix.
"""
function LinearAlgebra.lmul!(a::Number, K::AbstractKron)
    A = get_matrices(K)
    lmul!(a, A[1])
    return K
end

"""
    rmul!(K::KroneckerProduct, a::Number)

Scale a `KroneckerProduct` `K` inplace by a factor `a` by rescaling the
**right** matrix.
"""
function LinearAlgebra.rmul!(K::AbstractKron, a::Number)
    A = get_matrices(K)
    rmul!(A[end], a)
    return K
end

function Base.:*(a::Number, K::AbstractKron)
    A = get_matrices(K)
    KroneckerProduct(a * A[1], A[2:end]...)
end

function Base.:*(K::AbstractKron, a::Number)
    A = get_matrices(K)
    KroneckerProduct(A[1:end-1]..., A[end] * a)
end
