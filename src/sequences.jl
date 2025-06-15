using Base: Indices, @propagate_inbounds

export Interval, measure, SortedSequence, IncreasingSequence, NonDecreasingSequence
export IncreasingVector, NonDecreasingVector, IncreasingRange
export global_insert, insert!

"""
    SortedSequence{T}

Abstract type for sorted sequences of real numbers.
"""
abstract type SortedSequence{T<:Real} <: AbstractVector{T} end

"""
    IncreasingSequence{T}

Abstract type for increasing sequences of real numbers.
"""
abstract type IncreasingSequence{T<:Real} <: SortedSequence{T} end

"""
    Interval(a::T, b::T)

Create an interval from a point a to a point b.
"""
struct Interval{T<:Real} <: IncreasingSequence{T}
    a::T
    b::T
    function Interval(a, b)
        if !(a < b)
            throw(ArgumentError("Values need to be in increasing order."))
        end
        T = promote_type(typeof(a), typeof(b))
        return new{T}(T(a), T(b))
    end
end

Base.size(::Interval) = (2,)
Base.size(I::Interval, k) = size(I)[k]

function Base.getindex(I::Interval, i::Int)
    if i==1
        return I.a
    elseif i==2
        return I.b
    else
        throw(BoundsError("Index is out of range."))
    end
end

@inline measure(I::Interval) = I.b-I.a

@inline function (I::Interval)(x::T) where {T<:Real}
    return (1.0-x) * I.a + x * I.b
end

"""
    NonDecreasingSequence{T}

Abstract type for non-decreasing sequences of real numbers.
"""
abstract type NonDecreasingSequence{T<:Real} <: SortedSequence{T} end

Base.size(v::SortedSequence) = size(v.data)
Base.copy(v::SortedSequence) = deepcopy(v)

@inline @propagate_inbounds function Base.getindex(v::SortedSequence, i)
    v.data[i]
end


"""
    IncreasingRange{T<:Real}

Construct an increasing range of real numbers.
"""
struct IncreasingRange{T<:Real} <: IncreasingSequence{T}
    data::LinRange{T}
    function IncreasingRange(range::LinRange{T}) where T
        if !(range.start < range.stop)
            throw(ArgumentError("Values need to be in increasing order."))
        end
        return new{T}(range)
    end
end

function IncreasingRange(a, b, n::Int)
    return IncreasingRange(LinRange(a, b, n))
end

IncreasingRange(dom::Interval, n::Int) = IncreasingRange(dom.a, dom.b, n) 


IncreasingRange{T}(a, b, n::Int) where {T<:Real} = IncreasingRange(T(a), T(b), n)


"""
    IncreasingVector{T<:Real}

Construct a vector with an increasing set of real numbers.

# Examples:
The element type `T` is a subtype of `Real`. Hence it is
possible to make `IncreasingVector{Int64}` as well as `IncreasingVector{Float64}`
```jldoctest
julia> IncreasingVector([0,1,2,3])
4-element IncreasingVector{Int64}:
 0
 1
 2
 3

julia> IncreasingVector([0.0,1.0,2.0,3.0])
4-element IncreasingVector{Float64}:
 0.0
 1.0
 2.0
 3.0
```
It is also possible to extract the unique values of a NonDecreasingVector
into an IncreasingVector.
```jldoctest
julia> v = NonDecreasingVector([0.0,0.0,1.0,2.0,2.0,3.0,3.0])
7-element NonDecreasingVector{Float64}:
 0.0
 0.0
 1.0
 2.0
 2.0
 3.0
 3.0

julia> IncreasingVector(v)
4-element IncreasingVector{Float64}:
 0.0
 1.0
 2.0
 3.0
```
"""
struct IncreasingVector{T<:Real} <: IncreasingSequence{T}
    data::Vector{T}
    function IncreasingVector(v::Vector{T}, verify_ordering::Bool = true) where {T<:Real}
        if verify_ordering
            for i = 2:length(v)
                if !(v[i-1] < v[i])
                    throw(ArgumentError("Values need to be in increasing order."))
                end
            end
        end
        return new{T}(v)
    end
end

function IncreasingVector(x::IncreasingRange)
    return IncreasingVector(Vector(x.data), false)
end

IncreasingVector(x::IncreasingVector) = x
IncreasingVector(x) = IncreasingVector(Vector(x)) # fall-back method, inheriting from Vector

@inline @propagate_inbounds function Base.getindex(v::IncreasingVector, i::Int)
    v.data[i]
end

@inline @propagate_inbounds function Base.getindex(v::IncreasingVector, I)
    IncreasingVector(v.data[I], false)
end

"""
    NonDecreasingVector{T}

Construct a vector with a non-decreasing set of real numbers.

# Examples:
```jldoctest
julia> v = NonDecreasingVector([0.0,0.0,1.0,2.0,2.0])
5-element NonDecreasingVector{Float64}:
 0.0
 0.0
 1.0
 2.0
 2.0

julia> u, m = deconstruct_vector(v)
([0.0, 1.0, 2.0], [2, 1, 2], [1, 1, 2, 3, 3], [2, 3, 5])

julia> u
3-element IncreasingVector{Float64}:
 0.0
 1.0
 2.0

julia> NonDecreasingVector(u, m) == v
true
```
It is also possible to extract the unique values of a NonDecreasingVector
into an IncreasingVector.
```jldoctest
julia> v = NonDecreasingVector([0.0,0.0,1.0,2.0,2.0,3.0,3.0])
7-element NonDecreasingVector{Float64}:
 0.0
 0.0
 1.0
 2.0
 2.0
 3.0
 3.0

julia> IncreasingVector(v)
4-element IncreasingVector{Float64}:
 0.0
 1.0
 2.0
 3.0
```
"""
struct NonDecreasingVector{T<:Real} <: NonDecreasingSequence{T}
    data::Vector{T}
    function NonDecreasingVector(
        v::Vector{T},
        verify_ordering::Bool = true,
    ) where {T<:Real}
        if verify_ordering
            for i = 2:length(v)
                if (v[i-1] > v[i])
                    throw(ArgumentError("Values need to be in non-decreasing order."))
                end
            end
        end
        return new{T}(v)
    end
end

@inline @propagate_inbounds function Base.getindex(v::NonDecreasingVector, i::Int)
    v.data[i]
end

@inline @propagate_inbounds function Base.getindex(v::NonDecreasingVector, I)
    NonDecreasingVector(v.data[I], false)
end

function construct_vector(uvals::IncreasingSequence, mult::Vector{Int})
    u = construct_vector(uvals.data, mult)
    return NonDecreasingVector(u, false)
end

function deconstruct_vector(v::NonDecreasingVector)
    u, m, ia, ic = deconstruct_vector(v.data)
    return IncreasingVector(u, false), m, NonDecreasingVector(ia, false), IncreasingVector(ic, false)
end

function IncreasingVector(v::NonDecreasingSequence{T}) where {T}
    n = length(v)
    u = zeros(T, n)
    counter = 0
    for (index, item) in enumerate(Unique(v))
        u[index] = item[1]
        counter += 1
    end
    return IncreasingVector(u[1:counter], false)
end

function NonDecreasingVector(uvals::IncreasingSequence{T}, mult::Vector{Int}) where {T}
    return construct_vector(uvals, mult)
end

function NonDecreasingVector(uvals::Vector{T}, mult::Vector{Int}) where {T}
    u = construct_vector(uvals, mult)
    return NonDecreasingVector(u, true)
end

"""
    insert!(v, index, item)

The `Base.Insert!` function is extended to operate on `NonDecreasingVector`
and `IncreasingVector`. An argument-error will be thrown if the `item`
leads to an unsorted sequence.

# Examples:
```jldoctest
julia> v = NonDecreasingVector([1,2,3,4]);

julia> insert!(v, 3, 2)
5-element Vector{Int64}:
 1
 2
 2
 3
 4

julia> v = IncreasingVector([1,2,4,5]);

julia> insert!(v, 3, 3)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```
"""
function Base.insert!(
    v::NonDecreasingVector{T},
    index::Integer,
    item::T,
) where {T}
    if item > v[index] || (index > 1 && item < v[index-1])
        throw(ArgumentError("Values need to be in non-decreasing order."))
    end
    insert!(v.data, index, item)
end

function Base.insert!(v::IncreasingVector{T}, index::Integer, item::T) where {T}
    if !(item < v[index]) || (index > 1 && item <= v[index-1])
        throw(ArgumentError("Values need to be in increasing order."))
    end
    insert!(v.data, index, item)
end

"""
    global_insert(v, k)

Uniformly subdivide new values into an `IncreasingVector` or
`NonDecreasingVector`.

# Examples:
```jldoctest
julia> v = IncreasingVector([0.0,1.0,2.0])
3-element IncreasingVector{Float64}:
 0.0
 1.0
 2.0

julia> global_insert(v, 3)
9-element IncreasingVector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
 1.25
 1.5
 1.75
 2.0
```
When applied to a `NonDecreasingVector` only the non-zero length
intervals are subdivided
```jldoctest
julia> v = NonDecreasingVector([0.0,0.0,1.0,1.0,2.0])
5-element NonDecreasingVector{Float64}:
 0.0
 0.0
 1.0
 1.0
 2.0

julia> global_insert(v, 3)
11-element NonDecreasingVector{Float64}:
 0.0
 0.0
 0.25
 0.5
 0.75
 1.0
 1.0
 1.25
 1.5
 1.75
 2.0
```
"""
function global_insert(kts::NonDecreasingVector{T}, k::Int) where {T}

    # get unique knots and multiplicity
    ukts, umult, inda, indc = deconstruct_vector(kts)
    n = length(kts)
    m = length(ukts)

    # compute new knotvector
    newkts = zeros(T, n + k * (m - 1))
    j = 1
    for i = 1:m-1
        newkts[j:j+umult[i]-1] .= ukts[i]
        j += umult[i] - 1
        newkts[j:j+k+1] = LinRange(ukts[i], ukts[i+1], k + 2)
        j += k + 1
    end
    newkts[j:j+umult[m]-1] .= ukts[m]
    return NonDecreasingVector(newkts, false)
end

function global_insert(ukts::IncreasingVector{T}, k::Int) where {T}

    # get unique knots and multiplicity
    m = length(ukts)

    # compute new knotvector
    newkts = zeros(T, m + k * (m - 1))
    j = 1
    for i = 1:m-1
        newkts[j] = ukts[i]
        newkts[j:j+k+1] = LinRange(ukts[i], ukts[i+1], k + 2)
        j += k + 1
    end
    newkts[j] = ukts[m]
    return IncreasingVector(newkts, false)
end

function global_insert(v::IncreasingRange{T}, k::Int) where {T}
    range = v.data
    a = range.start
    b = range.stop
    n = range.len + range.lendiv * k
    IncreasingRange(a, b, n)
end
