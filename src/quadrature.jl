export QuadratureRule, standard_quadrature_rule, dimension_quadrature_rule, l2_error

IgaBase.standard_quadrature_rule(f, g::AbstractSpline) = IgaBase.QuadratureRule(TensorProduct(s -> PatchRule(s), g.space))

"""
    QuadratureRule(Q::TensorProduct{Dim,<:AbstractQuadrule{1}})

Generate a tensor product quadrature rule where the quadrature points are
stored as a CartesianProduct and the weights as a KroneckerProduct vector
allowing for easy evaluation and numerical integration based on Kronecker
products.
"""
function IgaBase.QuadratureRule(Q::TensorProduct{Dim,<:AbstractQuadrule{1}}) where Dim
    return IgaBase.QuadratureRule(CartesianProduct(q -> q.x[2:end-1], Q), KroneckerProduct(q -> q.w[2:end-1], Q, reverse=true))
end

IgaBase.dimension_quadrature_rule(x::CartesianProduct{Dim}, w) where {Dim} = Dim
