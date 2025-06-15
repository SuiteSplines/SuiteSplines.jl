export LinearIndex, MultiIndices, multiindex, linindex

import Combinatorics: binomial

"""
    Binomial(Val(n), Val(k)) = binomial(n,k)

`@generated` function for optimized computation of binomials.
"""
@generated Binomial(::Val{N}, ::Val{K}) where {N,K} = binomial(N,K)
Binomial(n::Integer, k::Integer) = Binomial(Val(n), Val(k))

"""
    LinearIndex(P::Degree, D::Dimension, I::Integer)

A linear index `I` of a simplicial discretization with Degree `P` and
dimension `D`.
"""
struct LinearIndex{P, D, I}
    LinearIndex(P::Degree, D::Dimension, I::Integer) = new{P,D,I}()
end
value(::LinearIndex{P,D,I}) where {P,D,I} = I

"""
    MultiIndex(::LinearIndex{P<:Degree, D<:Dimension, I<:Integer})

Return the multi-index associated with linear index `I` and simplicial
domain of degree `P` and dimension `D`.
"""
multiindex(::LinearIndex{P,2,I}) where {P,I} = (P+1-I, I-1) # step 1 in recursion 
@generated function multiindex(::LinearIndex{P,D,I}) where {P,D,I} # step n in recursion
    J = Binomial(Val(P+D-1),Val(D-1))
    Q = P+1
    while J>I-1
        J-=Binomial(Val(Q+D-3),Val(D-2))
        Q-=1
    end
    return (P-Q, (multiindex(LinearIndex(Q,D-1,I-J)))...)
end

"""
    linindex(i, j, k...)

Return the linear index of a multi-index label.
"""
@inline linindex(i::Integer, j::Integer) = j+1
@inline linindex(i::Integer, j::Integer, k::Integer) = Binomial(j+k+1, 2) + linindex(j,k)
@inline linindex(i::Integer, j::Integer, k::Integer, l::Integer) = Binomial(j+k+l+2,3) + linindex(j,k,l)
@inline linindex(i::Integer, j::Integer, k::Integer, l::Integer, m::Integer) = Binomial(j+k+l+m+3,4) + linindex(j,k,l,m)
@inline linindex(X::MultiIndex) = linindex(X...)

"""
    LinearIndex(::MultiIndex{X})

Compute the `LinearIndex` associated with a `MultiIndex`.
"""
@inline LinearIndex(X::MultiIndex) = LinearIndex(sum(X), length(X), linindex(X...))

"""
    MultiIndices(P::Degree, D::Dimension)
    MultiIndices(D::Dimension)
    MultiIndices(s::Simplex)
    MultiIndices(s::BezierSimplex)

Returns an iterator corresponding to a Bezier simplex of polynomial
degree `P` and dimension `D`.
"""
struct MultiIndices{P,D}
    n::Integer
end
MultiIndices(P::Degree, D::Dimension) = MultiIndices{P,D}(Binomial(Val(P+D-1), Val(D-1)))
MultiIndices(s::BezierSimplex{P,D}) where {P,D} = MultiIndices(P,D)
MultiIndices(D::Dimension) = MultiIndices(1, D)
MultiIndices(s::Simplex{D}) where {D} = MultiIndices(D)

Base.length(mi::MultiIndices) = mi.n

function Base.iterate(iter::MultiIndices{P,D}) where {P,D}
    return (multiindex(LinearIndex(P,D,1)), 1)
end

function Base.iterate(iter::MultiIndices{P,D}, state::Integer) where {P,D}
    state+=1
    if state < iter.n+1
        return (multiindex(LinearIndex(P,D,state)), state)
    end
end
