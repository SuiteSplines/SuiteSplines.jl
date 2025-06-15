using LinearAlgebra

using LinearAlgebra: checksquare

export det

box_exponent(m,n) = cld(m * n * (m-1) * (n-1), 4)
box_exponent(m,n,o...) = box_exponent(m,n) + box_exponent(m*n,o...)

"""
    det(K::BoxProduct)

Compute the determinant of a Box product.
"""
function LinearAlgebra.det(K::BoxProduct{M,T,2}) where {M,T}
    checksquare(K)
    if check(issquare, get_matrices(K))
        m = size(K,1)
        determinant = (-one(T))^box_exponent(get_row_sizes(K)...)
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
    inv(K::BoxProduct)

Compute the inverse of a Box product.
"""
function Base.inv(K::BoxProduct)
    checksquare(K)
    if check(issquare, get_matrices(K))
        return BoxProduct(inv, reverse(get_matrices(K)))
    else
        throw(SingularException(1))
    end
end

"""
    adjoint(K::BoxProduct)

Compute the adjoint of a Box product.
"""
function Base.adjoint(K::BoxProduct)
    return BoxProduct(adjoint, reverse(get_matrices(K)))
end

"""
    transpose(K::BoxProduct)

Compute the transpose of a Box product.
"""
function Base.transpose(K::BoxProduct)
    return BoxProduct(transpose, reverse(get_matrices(K)))
end

function Base.:*(A::BoxProduct{T,2}, B::BoxProduct{T,2}) where T
    if get_col_sizes(A) == reverse(get_row_sizes(B))
        # mixed product poperty
        return KroneckerProduct{2}(ntuple(k -> A.data[k] * B.data[end+1-k], order(A)))
    end
    @warn "Mismatch between A and D in (A ⊠ B)(C ⊠ D), using less efficient fall-back routines."
    return Matrix(A) *  Matrix(B)
end

function Base.:*(A::KroneckerProduct{T,2}, B::BoxProduct{T,2})  where T
    if get_col_sizes(A) == get_row_sizes(B)
        # mixed product poperty
        return BoxProduct{2}(ntuple(k -> A.data[k] * B.data[k], order(A)))
    end
    @warn "Mismatch between A and D in (A ⊗ B)(C ⊠ D), using less efficient fall-back routines."
    return Matrix(A) *  Matrix(B)
end

function Base.:*(A::BoxProduct{T,2}, B::KroneckerProduct{T,2}) where T
    if get_col_sizes(A) == reverse(get_row_sizes(B))
        # mixed product poperty
        return BoxProduct{2}(ntuple(k -> A.data[k] * B.data[end+1-k], order(A)))
    end
    @warn "Mismatch between A and D in (A ⊠ B)(C ⊗ D), using less efficient fall-back routines."
    return Matrix(A) *  Matrix(B)
end
