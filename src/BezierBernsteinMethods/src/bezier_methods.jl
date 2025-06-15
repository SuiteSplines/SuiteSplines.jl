export degree_elevation_operator, mass_matrix

"""
    degree_elevation_operator(p::Degree, d::Dimension)

Return a matrix that maps degrees of freedom of a Bezier polynomial of degree `p`
to the degrees of freedom corresponding to its order elevated Bezier polynomial.
"""
function degree_elevation_operator(p::Degree, d::Dimension, T::Type{<:Real} = Float64)
    m = dimsplinespace(p, d)       # spline space dim of Bezier patch
    n = dimsplinespace(p+1, d)     # spline space dim of higher order Bezier patch

    A = zeros(T,n,m)
    for (j,α) in enumerate(MultiIndices(p,d))
        for (k,δ) in enumerate(MultiIndices(d))
            i = linindex(α.+δ)
            A[i,j] = T(α[k]+1) / (p+1)
        end
    end
    return A
end

"""
    mass_matrix(S::AbstractSimplex{D}, p::Degree, q::Degree, T::Type{<:Real} = Type{Float64})

Compute Bernstein mass matrix for degree `(p,q)` and dimension `d`.
"""
function mass_matrix(S::AbstractSimplex{D}, p::Degree, q::Degree, T::Type{<:Real} = Float64) where D
    m = dimsplinespace(p,D)
    n = dimsplinespace(q,D)
    V = vol(S) * T(factorial(D-1)) / prod(p+q+1:p+q+D-1)
    M = zeros(T,n,m)

    for (j,α) in enumerate(MultiIndices(p,D))
        c = V * multinomial(α...)
        for (i,β) in enumerate(MultiIndices(q,D))
            M[i,j] = c * T(multinomial(β...)) / multinomial((α.+β)...)
        end
    end
    return M
end
