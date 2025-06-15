export @kronecker!, contract!, contract, mul!, mul, dot

import IgaBase: contract, contract!

macro kronecker!(ex)
    local op = esc(Val(ex.head))
    local E = esc(ex.args[1])
    local K = esc(ex.args[2].args[2])
    local X = esc(ex.args[2].args[3])
    return :(contract!($op, $E, $X, $K))
end

function contract(X::AbstractVector, V::AbstractVector)
    @tensor begin
        E = X[j] * V[j]
    end
    return E
end

function contract(X::AbstractVector{S}, K::KroneckerProduct{T,1,1}) where {S,T}
    B = get_matrices(K)
    @tensor begin
        E = X[j] * B[1][j]
    end
    return E
end

function contract(X::AbstractArray{S,2}, K::KroneckerProduct{T,1,2}) where {S,T}
    B = get_matrices(K)
    @tensor begin
        E = X[i,j] * B[2][i] * B[1][j]
    end
    return E
end

function contract(X::AbstractArray{S,3}, K::KroneckerProduct{T,1,3}) where {S,T}
    B = get_matrices(K)
    @tensor begin
        E = X[i,j,k] * B[3][i] * B[2][j] * B[1][k]
    end
    return E
end

function contract(X::AbstractArray{S,4}, K::KroneckerProduct{T,1,4}) where {S,T}
    B = get_matrices(K)
    @tensor begin
        E = X[i,j,k,l] * B[4][i] * B[3][j] * B[2][k] * B[1][l]
    end
    return E
end

function LinearAlgebra.dot(K::KroneckerProduct{T,1}, X::AbstractVector) where {T}
    m = reverse(get_row_sizes(K))
    return contract(reshape(X, m), K)
end

function LinearAlgebra.dot(X::AbstractVector, K::KroneckerProduct{T,1}) where {T}
    return dot(K,X)
end

# special cases: 1d (:=)
function contract!(::Val{:(=)}, E::AbstractVector, X::AbstractVector, K::KroneckerProduct{T,N,1}) where {T,N}
    E .= K * X
end

# special cases: 1d (:+=)
function contract!(::Val{:(+=)}, E::AbstractVector, X::AbstractVector, K::KroneckerProduct{T,N,1}) where {T,N}
    E .+= K * X
end

# special cases: 1d (:-=)
function contract!(::Val{:(-=)}, E::AbstractVector, X::AbstractVector, K::KroneckerProduct{T,N,1}) where {T,N}
    E .-= K * X
end

# special cases: 1d (:*=)
function contract!(::Val{:(*=)}, E::AbstractVector, X::AbstractVector, K::KroneckerProduct{T,N,1}) where {T,N}
    E .*= K * X
end

# special cases: 1d (:/=)
function contract!(::Val{:(/=)}, E::AbstractVector, X::AbstractVector, K::KroneckerProduct{T,N,1}) where {T,N}
    E ./= K * X
end

# general cases using TensorOperations
for N in 2:4
    local i = Any[ Symbol(:i_,i) for i in 1:N ]
    local j = Any[ Symbol(:j_,i) for i in 1:N ]
    local vars = Any[Expr(:ref, :(K.data[$(N+1-k)]), Symbol(:i,'_',k), Symbol(:j,'_',k)) for k in 1:N ]

    for OP in [:(=),:(+=),:(-=)]
        local ex = Expr(OP, Expr(:ref, :E, i...), Expr(:call, :*, Expr(:ref, :X, j...), vars...))
        local S = Val{OP}
        @eval function contract!(::$S, E::AbstractArray{T,$N}, X::AbstractArray{T,$N}, K::KroneckerProduct{T,2,$N}) where T
            @tensor begin
                $ex
            end
        end
    end
end

@inline function contract(op, X::AbstractArray{T,N}, K::KroneckerProduct{T,2,N}) where {T,N}
    sz = reverse(get_row_sizes(K))
    E = Array{promote_type(eltype(X),eltype(K))}(undef, sz...)
    return contract!(op, E, X, K)
end

function LinearAlgebra.mul!(y::AbstractVector, K::KroneckerProduct, x::AbstractVector)
    m, n = reverse(get_row_sizes(K)), reverse(get_col_sizes(K))
    contract!(Val(:(=)), reshape(y,m), reshape(x, n), K)
    return reshape(y,:)
end

function mul(K::KroneckerProduct, x::AbstractVector)
    m, n = reverse(get_row_sizes(K)), reverse(get_col_sizes(K))
    y = contract(Val(:(=)), reshape(x, n), K)
    return reshape(y, :)
end

Base.:*(K::KroneckerProduct{T,2,1}, x::AbstractVector) where T = K.data[1] * x
Base.:*(K::KroneckerProduct, x::AbstractVector) = mul(K,x)
