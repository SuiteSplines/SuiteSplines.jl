# define here any functions that should allow overloading inside
# the Feather eco-system

export Degree, Regularity, Dimension
export collect!, contract, contract!
export order, degree, orientation, direction, numbertype, getdata
export Even, Odd, EvenOrOdd, select_even_or_odd
export parameterless_typeof
export unzip, NotImplementedError
export find_singular_dimension, squeeze

"""
    Degree

Polynomial degree (typeallias of Integer).
"""
const Degree = Integer

"""
    Regularity

Regularity of the B-spline basis (typeallias of Integer).
"""
const Regularity = Integer

"""
    Dimension

Dimension (typeallias of Integer).
"""
const Dimension = Integer

function collect! end

function order end
function degree end
function orientation end
function direction end
function numbertype end
function getdata end
function contract end
function contract! end

abstract type EvenOrOdd  end
struct Even <: EvenOrOdd end
struct Odd  <: EvenOrOdd end

select_even_or_odd(m) = iseven(m) ? Even() : Odd()

unzip(c::NTuple{N,T}) where {N,T} = map(x -> getfield.(c, x), fieldnames(T))
unzip(c::Vector{T}) where {T} = map(x -> getfield.(c, x), fieldnames(T))

struct NotImplementedError <: Exception
    msg::AbstractString
end
NotImplementedError() = NotImplementedError("Not implemented.")

function Base.showerror(io::IO, err::NotImplementedError)
    print(io, "NotImplementedError: ")
    print(io, err.msg)
end

parameterless_typeof(::T) where T = T.name.wrapper

# find the first dimension that is singular
function find_singular_dimension(S::Int...)
    for k in 1:length(S)
        if S[k]==1
            return k
        end
    end
    return -1
end

# remove singular dimensions if there are any
function squeeze(A::AbstractArray)
    s = find_singular_dimension(size(A)...)
    if s ∈ 1:ndims(A)
        return dropdims(A, dims=s)
    end
    return A
end
