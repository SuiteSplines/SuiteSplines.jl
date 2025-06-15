"""
    TPEvaluationCache{Dim}

Type that caches spline basis functions as a `KroneckerProduct{Dim}` and a
grid of points as a `CartesianProduct{Dim}`. The type is mutable such that
the cached variables can be updated at runtime.
"""
mutable struct TPEvaluationCache{Dim} <: EvaluationCache{Dim}
    func::NTuple{Dim,Function}
    basis::KroneckerProduct
    grid::CartesianProduct{Dim}
    isinit::Bool
    function TPEvaluationCache(f::NTuple{Dim,Function}) where {Dim}
        eval = new{Dim}()
        eval.func = f
        eval.isinit = false
        return eval
    end
end

TPEvaluationCache(f...) = TPEvaluationCache(f)
TPEvaluationCache(f, iterable) = TPEvaluationCache((map(f, iterable)...,))
TPEvaluationCache(f, iterable...) = TPEvaluationCache((map(f, iterable...)...,))

# The cache is updated on initialization and afterwords when the gridpoints
# change
function IgaBase.update!(eval::TPEvaluationCache, grid::CartesianProduct)
    if isinitialized(eval) && (grid==eval.grid)
        return false # no update
    end
    eval.grid = grid
    eval.basis = KroneckerProduct((f,x) -> f(x), eval.func, eval.grid.data; reverse=true)
    eval.isinit = true
    return true
end
