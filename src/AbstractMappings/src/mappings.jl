export AbstractMapping, ScalarMapping, ScalarFunction, dimension, codimension, GeometricMapping, DifferentialForm
export types, boundary, ∂
export Gradient, Hessian, Normal, UnitNormal, Vol

import IgaBase: AbstractMapping, ScalarMapping, dimension, codimension, n_input_args, n_output_args, boundary

n_input_args(f::AbstractMapping) = dimension(f)

function n_input_args(first, args...)
    n = n_input_args(first)
    for f in args
        n_input_args(f) != n && throw(ArgumentError("Functors require the same number of input arguments."))
    end
    return n
end

function n_output_args(args...)
    return length(args)
end

function n_input_args(f::Function)::Int
    m = methods(f)
    n = length(m)
    @assert n == 1 "expected exactly one implemented method!"
    return first(m).nargs - 1
end

Base.size(f::AbstractMapping) = codimension(f)
Base.length(f::AbstractMapping) = prod(size(f))
Base.iterate(f::AbstractMapping) = iterate(f.data)
Base.iterate(f::AbstractMapping, state) = iterate(f.data, state)
Base.getindex(f::AbstractMapping, k...) = f.data[k...]

struct ScalarFunction{Dim} <: ScalarMapping{Dim}
    data::Function
    orientation::Int8
    function ScalarFunction(f::Function; orientation::Int=1)
        Dim = n_input_args(f)
        @assert f(zeros(Dim)...) isa Real "The provided function is non-scalar."
        return new{Dim}(f, orientation)
    end
end
@inline (f::ScalarFunction)(x...) = f.data(x...)

process_mapping_input(f) = f
process_mapping_input(f::Function) = ScalarFunction(f)

isa_domain(::Any) = false
isa_domain(::Interval) = true
function isa_domain(dom::CartesianProduct{Dim}) where Dim
    return dom.data isa NTuple{Dim,Interval}
end

IgaBase.domain(f::AbstractMapping) = f.domain

function check_mapping_arguments(args...) end

function types end

struct GeometricMapping{Dim,Codim,Dom,S<:Tuple} <: AbstractMapping{Dim,1,Codim}
    domain::Dom
    data::S
    orientation::Int8
    function GeometricMapping(domain, args...; orientation::Int=1)
        @assert isa_domain(domain) "First argument should be a domain type."

        # This does nothing right now:
        # neither NURBS nor TensorProductBsplines implement that
        # and the default (see this file, line 56) returns `nothing`
        check_mapping_arguments(args)
        
        Dim = n_input_args(args...)
        @assert Dim == ndims(domain)
        y = map(process_mapping_input, args)
        Codim = n_output_args(y...)
        Dom = typeof(domain)
        S = typeof(y)
        return new{Dim,Codim,Dom,S}(domain, y, orientation)
    end
end

function GeometricMapping(T::Type, args...; codimension::Int=1, orientation::Int=1)
    x = T(args...)
    y = ntuple(k -> similar(x), codimension-1)
    return GeometricMapping(domain(x), x, y...; orientation=orientation)
end

types(mapping::GeometricMapping) = typeof(mapping.data)

function get_property_imp(::Type, f, s) end

function Base.getproperty(f::GeometricMapping, s::Symbol)
    s === :domain && return getfield(f, :domain)
    s === :data && return getfield(f, :data)
    s === :orientation && return getfield(f, :orientation)
    return get_property_imp(types(f), f, s)
end

function IgaBase.boundary_imp(f::GeometricMapping, comp, dir)
    x = ntuple(k -> IgaBase.boundary_imp(f[k], f.domain, comp, dir), length(f))
    return GeometricMapping(IgaBase.boundary_imp(f.domain, comp, dir), x...; orientation=IgaBase.orientation(comp,dir))
end


IgaBase.boundary_imp(X::CartesianProduct{2}, comp, dir::Val{1}) = X.data[2]
IgaBase.boundary_imp(X::CartesianProduct{2}, comp, dir::Val{2}) = X.data[1]
IgaBase.boundary_imp(X::TensorProduct{2}, comp, dir::Val{1}) = X.data[2]
IgaBase.boundary_imp(X::TensorProduct{2}, comp, dir::Val{2}) = X.data[1]

IgaBase.boundary_imp(X::CartesianProduct, comp, dir) = CartesianProduct(IgaBase.boundary_imp(X.data, comp, dir)...)
IgaBase.boundary_imp(X::TensorProduct, comp, dir) = TensorProduct(IgaBase.boundary_imp(X.data, comp, dir)...)

# boundary implementation of ScalarFunction
for Comp in 1:2

    @eval function IgaBase.boundary_imp(f::ScalarFunction{2}, domain::CartesianProduct{2}, comp::Val{$Comp}, dir::Val{1})
        return ScalarFunction(y -> f.data(domain.data[1][$Comp], y); orientation=IgaBase.orientation(comp,dir))
    end
    @eval function IgaBase.boundary_imp(f::ScalarFunction{2}, domain::CartesianProduct{2}, comp::Val{$Comp}, dir::Val{2})
        return ScalarFunction(x -> f.data(x, domain.data[2][$Comp]); orientation=IgaBase.orientation(comp,dir))
    end

    # 3D
    @eval function IgaBase.boundary_imp(f::ScalarFunction{3}, domain::CartesianProduct{3}, comp::Val{$Comp}, dir::Val{1})
        return ScalarFunction((y,z) -> f.data(domain.data[1][$Comp],y,z); orientation=IgaBase.orientation(comp,dir))
    end
    @eval function IgaBase.boundary_imp(f::ScalarFunction{3}, domain::CartesianProduct{3}, comp::Val{$Comp}, dir::Val{2})
        return ScalarFunction((x,z) -> f.data(x,domain.data[2][$Comp],z); orientation=IgaBase.orientation(comp,dir))
    end
    @eval function IgaBase.boundary_imp(f::ScalarFunction{3}, domain::CartesianProduct{3}, comp::Val{$Comp}, dir::Val{3})
        return ScalarFunction((x,y) -> f.data(x,y,domain.data[3][$Comp]); orientation=IgaBase.orientation(comp,dir))
    end
end

struct Gradient{Dim,Codim,T<:AbstractMapping{Dim}} <: AbstractMapping{Dim,Dim,Codim}
    mapping::T
    function Gradient(f::AbstractMapping{Dim,1,1}) where {Dim}
        T = typeof(f)
        new{Dim,1,T}(f)
    end
    function Gradient(f::AbstractMapping{Dim,1,Codim}) where {Dim,Codim}
        T = typeof(f)
        new{Dim,Codim,T}(f)
    end
    function Gradient(f::AbstractMapping{Dim,Codim,1}) where {Dim,Codim}
        T = typeof(f)
        new{Dim,Codim,T}(f)
    end
end

struct Hessian{Dim,T<:AbstractMapping{Dim,1,1}} <: AbstractMapping{Dim,Dim,Dim}
    mapping::T
    function Hessian(f::AbstractMapping{Dim,1,1}) where Dim
        T = typeof(f)
        new{Dim,T}(f)
    end
end

struct UnitNormal{Dim,Codim,T<:GeometricMapping{Dim}} <: AbstractMapping{Dim,Codim,1}
    mapping::T
    function UnitNormal(f::GeometricMapping{Dim,Codim}) where {Dim,Codim}
        T = typeof(f)
        new{Dim,Codim,T}(f)
    end
end

struct Normal{Dim,Codim,T<:GeometricMapping{Dim}} <: AbstractMapping{Dim,Codim,1}
    mapping::T
    function Normal(f::GeometricMapping{Dim,Codim}) where {Dim,Codim}
        T = typeof(f)
        new{Dim,Codim,T}(f)
    end
end

struct Vol{Dim,T<:GeometricMapping{Dim}} <: AbstractMapping{Dim,1,1}
    mapping::T
    function Vol(f::GeometricMapping{Dim}) where {Dim}
        T = typeof(f)
        new{Dim,T}(f)
    end
end

@inline function Base.getindex(∇f::Gradient{1}, j::Int)
    return partial_derivative(∇f.mapping[j], 1)
end

@inline function Base.getindex(∇f::Gradient{Dim,1}, i::Int) where {Dim}
    return partial_derivative(∇f.mapping[1], ntuple(l -> Int(l==i), Dim))
end

@inline function Base.getindex(∇f::Gradient{Dim,Codim}, i::Int, j::Int) where {Dim,Codim}
    return partial_derivative(∇f.mapping[j], ntuple(l -> Int(l==i), Dim))
end

@inline function Base.getindex(Δf::Hessian{Dim}, i::Int, j::Int, k...) where Dim
    return partial_derivative(Δf.mapping[1], ntuple(l -> Int(l==i)+Int(l==j), Dim))
end