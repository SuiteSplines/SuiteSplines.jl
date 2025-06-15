using LinearAlgebra: Algorithm, default_svd_alg, checksquare

export eigen, svd, svd!, cholesky, cholesky!

##region KroneckerProduct factorizations
function LinearAlgebra.:\(K::KroneckerProduct{T,2}, c::AbstractVector) where T
    size(K, 1) != length(c) && throw(DimensionMismatch("size(K, 1) != length(c)"))
    return inv(K) * c
end

"""
    eigen(K::KroneckerProduct)

If the matrices of a `KroneckerProduct` instance are square, performs the eigenvalue
decompositon on them and return an `Eigen` type.
"""
function LinearAlgebra.eigen(K::KroneckerProduct)
    map(checksquare, get_matrices(K))
    λ, ϕ = unzip(map(eigen, get_matrices(K)))
    return Eigen(KroneckerProduct(λ...), KroneckerProduct(ϕ...))
end

function LinearAlgebra.svd(K::KroneckerProduct{T}; full::Bool = false, alg::Algorithm=default_svd_alg(K)) where T
    u, s, vt = unzip(map(A -> svd(A, full=full, alg=alg), get_matrices(K)))
    U, S, Vt = KroneckerProduct(u...), KroneckerProduct(s...), KroneckerProduct(vt...)
    return SVD{T, T, typeof(U)}(U, S, Vt)
end

function LinearAlgebra.svd!(K::KroneckerProduct; full::Bool = false, alg::Algorithm=default_svd_alg(K))
    U, S, Vt = unzip(map(A -> svd!(A, full=full, alg=alg), get_matrices(K)))
    return SVD(KroneckerProduct(U...), KroneckerProduct(S...),  KroneckerProduct(Vt...))
end

function LinearAlgebra.cholesky(A::KroneckerProduct; check = true)
    factors, uplo, info = unzip(map(a -> cholesky(a; check = check), get_matrices(A)))
    A = KroneckerProduct{2}(factors)
    Cholesky{eltype(A),typeof(A)}(A, uplo[1], info[1])
    # return Cholesky(KroneckerProduct{2}(factors), uplo[1], info[1])
end

function LinearAlgebra.cholesky!(A::KroneckerProduct; check = true)
    factors, uplo, info = unzip(map(a -> cholesky!(a; check = check), get_matrices(A)))
    return Cholesky(KroneckerProduct{2}(factors), uplo[1], info[1])
end

function Base.getproperty(C::Cholesky{T,<:KroneckerProduct{T, 2}}, d::Symbol) where T
    Cfactors = get_matrices(getfield(C, :factors))
    Cuplo    = getfield(C, :uplo)
    info     = getfield(C, :info)
    if d === :U
        if Cuplo === LinearAlgebra.char_uplo(d)
            return KroneckerProduct{2}(ntuple(k -> UpperTriangular(Cfactors[k]), length(Cfactors)))
        else
            return KroneckerProduct{2}(ntuple(k -> UpperTriangular(copy(Cfactors[k]')), length(Cfactors)))
        end
    elseif d === :L
        if Cuplo === LinearAlgebra.char_uplo(d)
            return KroneckerProduct{2}(ntuple(k -> LowerTriangular(Cfactors[k]), length(Cfactors)))
        else
            return KroneckerProduct{2}(ntuple(k -> LowerTriangular(copy(Cfactors[k]')), length(Cfactors)))
        end
    elseif d === :UL
        if Cuplo === 'U'
            return KroneckerProduct{2}(ntuple(k -> UpperTriangular(Cfactors[k]), length(Cfactors)))
        else
            return KroneckerProduct{2}(ntuple(k -> LowerTriangular(Cfactors[k]), length(Cfactors)))
        end
        return getfield(C, d)
    end
    getfield(C, d)
end




# function LinearAlgebra.lu(K::KroneckerProduct, pivot::Union{Val{false}, Val{true}} = Val(true); check = true)
#     factors, ipiv, info = unzip(map(A -> lu(A, pivot; check = check), get_matrices(K)))
#     n = map(length, ipiv)
#     return LU(KroneckerProduct(factors...), reshape(1:prod(n), n[end:-1:1])[ipiv[end:-1:1]...][:], info[1])
# end
#
# function LinearAlgebra.lu!(K::KroneckerProduct, pivot::Union{Val{false}, Val{true}} = Val(true); check = true)
#     factors, ipiv, info = unzip(map(A -> lu!(A, pivot; check = check), get_matrices(K)))
#     n = map(length, ipiv)
#     return LU(KroneckerProduct(factors...), reshape(1:prod(n), n[end:-1:1])[ipiv[end:-1:1]...][:], info[1])
# end
#
# function Base.getproperty(F::LU{T,KroneckerProduct{M, T, 2}}, d::Symbol) where {T,M}
#     m, n = size(F)
#     if d === :L
#         L = map(A -> LowerTriangular(copy(A)), get_matrices(F.factors))
#         for k in 1:M
#             for i in 1:size(L[k],1)
#                 L[k][i,i] = one(T)
#             end
#         end
#         return KroneckerProduct(L...)
#     elseif d === :U
#         U = map(A -> UpperTriangular(A), get_matrices(F.factors))
#         return KroneckerProduct(U...)
#     elseif d === :p
#         return LinearAlgebra.ipiv2perm(getfield(F, :ipiv), m)
#     elseif d === :P
#         return Matrix{T}(I, m, m)[:,invperm(F.p)]
#     else
#         getfield(F, d)
#     end
# end

##endregion


##region KroneckerSum factorizations

function LinearAlgebra.eigen(K::KroneckerSum)
    λ, ϕ = unzip(map(eigen, get_matrices(K)))
    return Eigen(KroneckerSum(λ...), KroneckerProduct(ϕ...))
end


##endregion
