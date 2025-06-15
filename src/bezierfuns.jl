export dimsplinespace, grevillepoint, bernsteinfuns, bernstein, decasteljau!

import Combinatorics: multinomial

"""
    dimsplinespace(p::Degree, d::Dimension)

Compute the dimension of the Bezier spline space on a simplex of
degree `p` and numer of vertices `d`.
"""
@inline dimsplinespace(p::Degree, d::Dimension) = Binomial(p+d-1, d-1)

"""
    grevillepoint(X::MultiIndex)

Compute Greville point associated with multi-index label `X`.
"""
@inline grevillepoint(X::MultiIndex) = SVector(broadcast(/, X, sum(X)))

"""
    bernstein(X::MultiIndex{D}, λ::BPoint{D,T})

Compute single Bernstein basis function corresponding to multi-index `X` at
barycentric point `λ`.
"""
@generated function bernsteinfuns(::Val{X}, λ::BPoint{D,T}) where {X,D,T}
    @assert length(X)==D
    tmp = T(multinomial(X...))
    quote
        s = $tmp
        for (k,x) in enumerate(X)
            s *= λ[k]^x
        end
        return s
    end
end
bernsteinfuns(X::MultiIndex{D}, λ::BPoint{D,T}) where {D,T} = bernsteinfuns(Val(X), λ)

"""
    bernstein(p::Degree, control::AbstractVector{T}, λ::BPoint{D,T})

Compute Bezier polynomial at barycentric point `λ` with the degrees of freedom in `control`.
"""
function bernstein(p::Degree, control::AbstractVector{T}, λ::BPoint{D,T}) where {D,T}
    @assert length(control)==dimsplinespace(p,D)
    s = T(0)
    for (i,μ) in enumerate(MultiIndices(p,D))
        s += control[i] * bernsteinfuns(μ, λ)
    end
    return s
end

"""
    step_decasteljau!(control::Vector{T}, λ::BPoint{D,T}, MultiIndices)

Computes a single step of the recursive DeCasteljau algorithm on a `(D-1)`-dimensional
simplex.
"""
function step_decasteljau!(control::Vector{T}, λ::BPoint{D,T}, multiindices) where {D,T}
    for μ in multiindices
        s = T(0)
        for (j,δ) in enumerate(MultiIndices(D))
            s += λ[j] * control[value(LinearIndex(μ .+ δ))]
        end
        control[value(LinearIndex(μ))] = s;
    end
end

"""
    decasteljau!(control::Vector{T}, λ::BPoint{D,T}, MultiIndices)

Classic DeCasteljau algorithm (can compute derivatives too!) on a `(D-1)`-dimensional
simplex.
"""
function decasteljau!(p::Degree, control::AbstractVector{T}, λ::BPoint{D,T}) where {D,T}
    @assert length(control)==dimsplinespace(p,D)
    for r in 1:p
        step_decasteljau!(control, λ, MultiIndices(p-r, D))
    end
    return control[1]
end

function decasteljau!(p::Degree, control::AbstractVector{T}, λ::BPoint{D,T}, δλ::BVector{D,T}) where {D,T}

    # Loop over barycentric directions
    k = length(δλ)
    c = T(1)
    for r in 1:k
        step_decasteljau!(control, δλ[r], MultiIndices(p-r, D))
        c *= p + 1 - r;
    end

    # Loop over barycentric coordinates
    decasteljau!(p-k, control, λ)

    return c * control[1]
end
