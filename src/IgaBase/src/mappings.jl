export AbstractMapping, ScalarMapping, AbstractSpline, dimension, codimension, domain, n_input_args, n_output_args
export EvaluationCache, isinitialized, update!
export differentiate, partial_derivative
export restrict, restrict_to, @restrict
export dir, comp, boundary, ∂, Boundary, num_boundary_components, boundary_element_type, restriction

"""
    AbstractMapping{Dim,S1,S2}

Abstract type representing a mapping with an `Dim` dimensional
domain, taking values in an `S1 x S2` codimensional space.
"""
abstract type AbstractMapping{Dim,S1,S2} end

"""
    ScalarMapping{Dim} <: AbstractMapping{Dim,1,1}

Abstract subtype of `AbstractMapping`` that represents scalar-valued
mappings.
"""
abstract type ScalarMapping{Dim} <: AbstractMapping{Dim,1,1} end

"""
    AbstractSpline{Dim} <: ScalarMapping{Dim}

Abstract subtype of `ScalarMapping` that represents a AbstractSpline
function.
"""
abstract type AbstractSpline{Dim} <: ScalarMapping{Dim} end


dimension(::AbstractMapping{Dim}) where {Dim} = Dim
codimension(::AbstractMapping{Dim,S1,S2}) where {Dim,S1,S2} = (S1,S2)
codimension(f::AbstractMapping, k::Int) = codimension(f)[k]

function domain end
function n_input_args end
function n_output_args end
function differentiate end
function partial_derivative end

abstract type EvaluationCache{Dim} end

# Concrete subtypes of EvaluationCache have a mutable struct layout that
# looks like
#
# mutable struct MyEvaluationCache <: EvaluationCache{Dim}
#     func::F
#     basis::S
#     grid::T
#     isinit::Bool
#     function EvaluationCache(f::Function)
#         eval = new()
#         eval.func = f
#         eval.isinit = false
#         return eval
#     end
# end

isinitialized(eval::EvaluationCache) = eval.isinit
Base.ndims(::EvaluationCache{Dim}) where Dim = Dim

# A minimal implementation implements the following functionality
update!(eval::EvaluationCache, x) = throw(NotImplementedError("update! is not implemented."))

# implement restrict to get a presonalized view on your data
function restrict end

function restrict_imp end

# this piece of code is not imported from Base such that it works with
# all versions of Julia
include("replace_ref_begin_end.jl")

# macro that provides a view of an AbstractArray, e.g.
# @restrict A[:,:,1]
# requires implementation of ``restrict(A, I)``
macro restrict(ex)
    if Meta.isexpr(ex, :ref)
        ex = custom_replace_ref_begin_end!(ex)
        if Meta.isexpr(ex, :ref)
            ex = Expr(:call, restrict, ex.args...)
        else # ex replaced by let ...; foo[...]; end
            @assert Meta.isexpr(ex, :let) && Meta.isexpr(ex.args[2], :ref)
            ex.args[2] = Expr(:call, restrict, ex.args[2].args...)
        end
        Expr(:&&, true, esc(ex))
    else
        throw(ArgumentError("Invalid use of @restrict macro: argument must be a reference expression A[...]."))
    end
end

# similar as restrict, but provides a view on the boundary
# data of an array
restrict_to(A::AbstractArray; side) = restrict_imp(A, Val(comp(side)), Val(dir(side)))

# restrict_to implementation of AbstractArray
for Comp in 1:2

    index = Comp==1 ? :1 : :end

    # 2D case
    @eval restrict_imp(A::AbstractArray{T,2}, ::Val{$Comp}, ::Val{1}) where T = @restrict A[$index, :]
    @eval restrict_imp(A::AbstractArray{T,2}, ::Val{$Comp}, ::Val{2}) where T = @restrict A[:, $index]

    # 3D case
    @eval restrict_imp(A::AbstractArray{T,3}, ::Val{$Comp}, ::Val{1}) where T = @restrict A[$index, :, :]
    @eval restrict_imp(A::AbstractArray{T,3}, ::Val{$Comp}, ::Val{2}) where T = @restrict A[:, $index, :]
    @eval restrict_imp(A::AbstractArray{T,3}, ::Val{$Comp}, ::Val{3}) where T = @restrict A[:, :, $index]
end


# functionality for boundaries
function boundary end
function boundary_imp end
function num_boundary_components end
function boundary_element_type end

"""
    boundary(s; direction=1, component=1)
    ∂(s; direction=1, component=1)

Return boundary `component` in `direction` of the object `s`.
"""
boundary(s, component::Int, direction::Int) = IgaBase.boundary_imp(s, Val(component), Val(direction))

@inline IgaBase.boundary_imp(s, comp, dir) = IgaBase.boundary_imp(s, domain(s), comp, dir) 

@inline dir(k::Int)  = div(k+1, 2)      # direction given linear index
@inline comp(k::Int) = mod(k+1, 2)+1    # component given linear index

function boundary(s, k::Int)
    return boundary(s, comp(k), dir(k))
end

∂(args...) = boundary(args...)

# compute the boundary of datastructures defined on a rectangular domain.

# 2D
IgaBase.boundary_imp(data::NTuple{2}, comp, ::Val{1}) = (data[2],)
IgaBase.boundary_imp(data::NTuple{2}, comp, ::Val{2}) = (data[1],)

# 3D
IgaBase.boundary_imp(data::NTuple{3}, comp, ::Val{1}) = (data[2], data[3])
IgaBase.boundary_imp(data::NTuple{3}, comp, ::Val{2}) = (data[1], data[3])
IgaBase.boundary_imp(data::NTuple{3}, comp, ::Val{3}) = (data[1], data[2])

# Boundary implementation of AbstractArray
for Comp in 1:2

    index = Comp==1 ? :1 : :end

    # 2D case
    @eval boundary_imp_array(data::AbstractArray{T,2}, ::Val{$Comp}, ::Val{1}) where T = @view data[$index, :]
    @eval boundary_imp_array(data::AbstractArray{T,2}, ::Val{$Comp}, ::Val{2}) where T = @view data[:, $index]

    # 3D case
    @eval boundary_imp_array(data::AbstractArray{T,3}, ::Val{$Comp}, ::Val{1}) where T = @view data[$index, :, :]
    @eval boundary_imp_array(data::AbstractArray{T,3}, ::Val{$Comp}, ::Val{2}) where T = @view data[:, $index, :]
    @eval boundary_imp_array(data::AbstractArray{T,3}, ::Val{$Comp}, ::Val{3}) where T = @view data[:, :, $index]
end

@inline IgaBase.boundary_imp(data::AbstractArray, comp, dir) = boundary_imp_array(data, comp, dir)

# similar to boundary
restriction(data::AbstractArray, side) = boundary_imp_array(data, Val(comp(side)), Val(dir(side)))

# compute the orientation associated with taking the boundary of 
# datastructures defined on a rectangular domain.
for Comp in 1:2
    for Dir in 1:3
        @eval IgaBase.orientation(::Val{$Comp}, ::Val{$Dir}) = (-1)^($Comp+$Dir-1)
    end
end

"""
    Boundary{S<:AbstractMapping}

Iterator over the boundary components of a mapping of type `S`.
"""
struct Boundary{S<:AbstractMapping}
    data::S
end

# the following functions need to be implemented in order to use the
# Boundary iterator
Base.length(b::Boundary) = num_boundary_components(b.data)
Base.eltype(b::Boundary) = boundary_element_type(b.data)

function Base.iterate(iter::Boundary)
    return boundary(iter.data, 1), 1
end

function Base.iterate(iter::Boundary, state)
    if state < length(iter)
        state+=1
        return boundary(iter.data, state), state
    end
end