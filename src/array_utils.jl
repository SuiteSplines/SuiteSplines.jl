export construct_vector, deconstruct_vector, Unique

"""
    construct_vector(u, m)

Construct a vector by repeating the elements in `u` `m` times.

# Examples:
```jldoctest
julia> u = [0.0,1.0,2.5,3.0]
4-element Vector{Float64}:
 0.0
 1.0
 2.5
 3.0

julia> m = [2,1,2,3]
4-element Vector{Int64}:
 2
 1
 2
 3

julia> construct_vector(u, m)
8-element Vector{Float64}:
 0.0
 0.0
 1.0
 2.5
 2.5
 3.0
 3.0
 3.0
```
"""
function construct_vector(uvals::AbstractVector{T}, mult::Vector{Int}) where {T}
    @assert length(uvals) == length(mult) "Sizes of input argumens are not consistent."
    n = sum(mult)
    v = similar(uvals, n)
    ind = 1
    for i = 1:length(uvals)
        v[ind:ind+mult[i]-1] = fill(uvals[i], mult[i])
        ind = ind + mult[i]
    end
    return v
end

# count the number of times that v[index] is repeated, starting
# from index and counting in the direction of increasing indices
function count_multiplicity_down(v::AbstractVector, index)
    n = length(v)
    save = v[index]
    k = index + 1
    while !(k > n) && (v[k] == save)
        k += 1
    end
    return k - index
end

# count the number of times that v[index] is repeated, starting
# from index and counting in the direction of decreasing indices
function count_multiplicity_up(v::AbstractVector, index)
    n = length(v)
    save = v[index]
    k = index - 1
    while (k > 0) && (v[k] == save)
        k -= 1
    end
    return index - k
end

"""
    Unique(v)

Convenience iterator that lazily returns a tuple with the unique
consecutive values and their multiplicities.

# Examples:
```jldoctest
julia> v = [1,1,3,4,5,5,5,6,6,7];

julia> for item in Unique(v)
           @show item
       end
item = (1, 2)
item = (3, 1)
item = (4, 1)
item = (5, 3)
item = (6, 2)
item = (7, 1)
```
"""
struct Unique{T}
    data::AbstractVector{T}
end

Base.eltype(v::Unique) = Tuple{eltype(v.data),Integer}

function Base.iterate(iter::Unique)
    v = iter.data
    value = v[1]
    mult = count_multiplicity_down(v, 1)
    return ((value, mult), mult + 1)
end

function Base.iterate(iter::Unique, state::Integer)
    v = iter.data
    if !(state > length(v))
        value = v[state]
        mult = count_multiplicity_down(v, state)
        state += mult
        return ((value, mult), state)
    end
end

function Base.length(v::Unique)
    count = 0
    for item in v
        count += 1
    end
    return count
end

"""
    deconstruct_vector(v)

Decompose a sequence into a new sequence, a multiplicity vector,
and vectors specifying indexing into the sequence.

# Examples:
```jldoctest
julia> v = [0.0,0.0,1.0,2.0,2.0,3.0,3.0,3.0];

julia> u, m, ia, ic = deconstruct_vector(v)
([0.0, 1.0, 2.0, 3.0], [2, 1, 2, 3], [1, 1, 2, 3, 3, 4, 4, 4], [2, 3, 5, 8])

julia> construct_vector(u,m)==v
true

julia> u[ia]==v
true

julia> v[ic]==u
true

```
"""
function deconstruct_vector(v::AbstractVector{T}) where {T}

    n = length(v)
    mult = ones(Int64, n)
    vals = zeros(T, n)
    ia = zeros(Int64, n)
    ic = zeros(Int64, n)

    j1 = j2 = 0
    for (i, (u,μ)) in enumerate(Unique(v))

        # unique sequential values and multiplicity
        vals[i] = u
        mult[i] = μ

        # indexing vectors
        j2 = j1 + μ
        ia[j1+1:j2] .= i
        ic[i] = j2

        # initialize next iteration
        j1 = j2
    end
    last = ia[j2]

    return vals[1:last], mult[1:last], ia, ic[1:last]
end
