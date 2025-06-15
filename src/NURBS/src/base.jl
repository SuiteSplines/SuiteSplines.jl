export Nurbs, boundary, domain
export transform_to_projective_space!, transform_from_projective_space!

struct Nurbs{Dim, S<:AbstractSpline{Dim}} <: AbstractSpline{Dim}
    splinefun::S
    weightfun::S
    function Nurbs(splinefun::S, weightfun::S) where S
        @assert dimension(splinefun) == dimension(weightfun)
        @assert sum(splinefun.ders) === sum(weightfun.ders) == 0
        Dim = dimension(splinefun)
        return new{Dim,S}(splinefun, weightfun)
    end
end

function Base.getproperty(f::Nurbs, s::Symbol)
    s === :space && return getfield(getfield(f, :splinefun), :space)
    s === :coeffs && return getfield(getfield(f, :splinefun), :coeffs)
    s === :weights && return getfield(getfield(f, :weightfun), :coeffs)
    s === :ders && return getfield(getfield(f, :splinefun), :ders)
    s === :cache && return getfield(getfield(f, :splinefun), :cache)
    s === :orientation && return getfield(getfield(f, :splinefun), :orientation)
    return getfield(f, s)
end

function Base.propertynames(f::Nurbs)
    return (:space, :coeffs, :weights, :cache, :orientation)
end

UnivariateSplines.controlpoints(nurbs::Nurbs) = nurbs.coeffs ./ nurbsweights(nurbs)

function nurbsweights(nurbs::Nurbs{1})
    V = nurbs.weightfun.space
    W = nurbs.splinefun.space

    V == W && return nurbs.weights

    C = UnivariateSplines.two_scale_operator(V.p, V.U, W.p, W.U)
    weights = C * nurbs.weights
end

function nurbsweights(nurbs::Nurbs)
    V = nurbs.weightfun.space
    W = nurbs.splinefun.space

    V == W && return nurbs.weights

    f = (V, W) -> UnivariateSplines.two_scale_operator(V.p, V.U, W.p, W.U)
    C = KroneckerProduct(f, V, W; reverse=true)
    weights = zeros(size(W)...)
    @kronecker! weights = C * nurbs.weights
end

function Base.show(io::IO, spline::Nurbs{1})
    S = spline.space
    println(io, "Nurbs curve defined on splinespace")
    println(io, "$S")
end

function Base.show(io::IO, spline::Nurbs{Dim}) where {Dim}
    S = spline.space
    println(io, "Nurbs{$Dim} defined on splinespaces")
    for s in S
        println(io, "$s")
    end
end

function Nurbs(space, coeffs=zeros(numbertype(space),size(space)), weights=ones(numbertype(space),size(space)))
    S = spline_type(space)
    splinefun = S(space, coeffs)
    weightfun = S(space, weights, splinefun.cache; ders=splinefun.ders)
    return Nurbs(splinefun, weightfun)
end

function Base.similar(spline::Nurbs; coeffs=similar(spline.coeffs))
    return Nurbs(similar(spline.splinefun; coeffs=coeffs), spline.weightfun)
end

spline_type(::SplineSpace) = Bspline
spline_type(::TensorProduct{Dim,<:SplineSpace}) where Dim = TensorProductBspline

IgaBase.domain(spline::Nurbs) = IgaBase.domain(spline.space)

function transform_to_projective_space!(nurbs::Nurbs)
    nurbs.coeffs .*= nurbs.weights
end

function transform_from_projective_space!(nurbs::Nurbs)
    nurbs.coeffs ./= nurbsweights(nurbs)
end

# overloading of AbstractMappings.Evaluation.evalkernel! which enables computation
# on CartesianProduct grids of points using the @evaluation macro.
function AbstractMappings.evalkernel_imp!(op::Val{:(=)}, y, x, R::Nurbs)
    AbstractMappings.evalkernel_imp!(op, y, x, R.splinefun)
    w = similar(y)
    AbstractMappings.evalkernel_imp!(Val(:(=)), w, x, R.weightfun)
    _eval_nurbs!(y,w)
    return y
end

function IgaBase.boundary_imp(nurbs::Nurbs, domain, component, direction)
    s = IgaBase.boundary_imp(nurbs.splinefun, domain, component, direction)
    w = IgaBase.boundary_imp(nurbs.weightfun, domain, component, direction)
    return Nurbs(s, w)
end

# evaluate nurbs evaluation inplace
function _eval_nurbs!(Y::AbstractArray{T,N}, W::AbstractArray{T,N}) where {T,N}
    Y ./= W
end
