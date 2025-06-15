export AbstractQuadrule, QuadratureRule, GaussianRule, Legendre, Lobatto
export dimension_quadrature_rule, standard_quadrature_rule

abstract type AbstractQuadrule{Dim} end

Base.length(Q::AbstractQuadrule) = length(Q.x)
Base.size(Q::AbstractQuadrule) = size(Q.x)
Base.size(Q::AbstractQuadrule, k) = size(Q.x, k)
Base.copy(Q::AbstractQuadrule) = Base.deepcopy(Q)

# back-up method
dimension_quadrature_rule(x::AbstractVector, w) = 1
dimension_quadrature_rule(x::NTuple{Dim}, w) where {Dim} = Dim


function standard_quadrature_rule end

struct QuadratureRule{Dim,X,W} <: AbstractQuadrule{Dim}
    x::X
    w::W
    function QuadratureRule(x::X, w::W) where {X,W}
        Dim = dimension_quadrature_rule(x, w)
        return new{Dim,X,W}(x, w)
    end
end

# traits
abstract type GaussianRule end
struct Legendre <: GaussianRule end
struct Lobatto <: GaussianRule end