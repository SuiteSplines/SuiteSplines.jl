export BoxProduct, box, ⊠

"""
    BoxProduct{M,T,N} <: AbstractKron{M,T,N}

The box product is defined in the following paper.

Olsen, Peder A., Steven J. Rennie, and Vaibhava Goel. "Efficient automatic
differentiation of matrix functions." In Recent Advances in Algorithmic
Differentiation, pp. 71-81. Springer, Berlin, Heidelberg, 2012.
"""
struct BoxProduct{T,N,M,S<:Tuple} <: AbstractKron{T,N}
    data::S
    sz::Tuple{Int,Int}
    function BoxProduct{N}(A::S) where {N, S<:Tuple}
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

BoxProduct(A::AbstractVecOrMat...; reverse=false) = BoxProduct{2}((reverse==false) ? A : Base.reverse(A))
BoxProduct(A::AbstractVector...; reverse=false) = BoxProduct{1}((reverse==false) ? A : Base.reverse(A))
BoxProduct(f, iterable; reverse=false) = BoxProduct(map(f, iterable)...; reverse=reverse)
BoxProduct(f, iterable...; reverse=false) = BoxProduct(map(f, iterable...)...; reverse=reverse)

Base.size(K::BoxProduct) = K.sz
Base.size(K::BoxProduct, i) = K.sz[i]
Base.size(K::BoxProduct{T,1}) where T = (K.sz[1],)

⊠(A::AbstractVecOrMat, B::AbstractVecOrMat) = BoxProduct(A, B)
⊠(A::BoxProduct, B::AbstractVecOrMat) = BoxProduct(get_matrices(A)...,B)
⊠(A::AbstractVecOrMat, B::BoxProduct) = BoxProduct(A, get_matrices(B)...)

function Base.getindex(K::BoxProduct{T,2,M}, i1::Integer, i2::Integer) where {T,M}
    ia, ib = zeros(Int,M), zeros(Int,M)
    value = T(1)
    m, n = size(K)
    for k in 1:M
        l = M+1-k
        ia[k], i1, m = get_kronecker_indices(i1, m, size(K.data[k],1))
        ib[l], i2, n = get_kronecker_indices(i2, n, size(K.data[l],2))
    end
    for k in 1:M
        value *= K.data[k][ia[k], ib[k]]
    end
    return value
end

"""
    box(A, B)
    box(A, B...)

Compute the box product, which is defined as a perfect shuffle of a
kronecker product.
"""
function box end

function box(A)
    return Matrix(A)
end

function box(A, B)
    M, N = map(prod, unzip((size(A), size(B))))
    K = Matrix{promote_type(eltype(A),eltype(B))}(undef, M, N)
    for J in 1:N
        for I in 1:M
            m, n = M, N
            i1, i2 = I, J
            ia, i1, m = get_kronecker_indices(i1, m, size(A,1)); i = ia
            ia, i1, m = get_kronecker_indices(i1, m, size(B,1)); k = ia
            ib, i2, n = get_kronecker_indices(i2, n, size(B,2)); j = ib
            ib, i2, n = get_kronecker_indices(i2, n, size(A,2)); l = ib
            K[I,J] = A[i, l] * B[k, j]
        end
    end
    return K
end

box(A, B, C...) = box(box(A,B), C...)

Base.collect(K::BoxProduct) = box(get_matrices(K)...)
