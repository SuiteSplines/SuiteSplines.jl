export AbstractRefinement, hRefinement, pRefinement, kRefinement, hpRefinement
export refine

abstract type AbstractRefinement end

struct hRefinement <: AbstractRefinement
    h::Int
    function hRefinement(h)
        @assert h>0 "k should be an integer greater than zero."
        return new(h)
    end
end

struct pRefinement <: AbstractRefinement
    p::Int
    function pRefinement(p)
        @assert p>0 "Δp should be an integer greater than zero."
        return new(p)
    end
end

struct kRefinement <: AbstractRefinement
    h::Int
    p::Int
    function kRefinement(h, p)
        @assert p>0 "Δp should be an integer greater than zero."
        @assert h>0 "k should be an integer greater than zero."
        return new(h, p)
    end
end

struct hpRefinement <: AbstractRefinement
    h::Int
    p::Int
    function hpRefinement(h, p)
        @assert p>0 "Δp should be an integer greater than zero."
        @assert h>0 "k should be an integer greater than zero."
        return new(h, p)
    end
end

function refine(mapping::ScalarMapping; method::AbstractRefinement)
    return refine_imp(mapping, method)
end

function refine(mapping::AbstractMapping; method::AbstractRefinement)
    T = parameterless_typeof(mapping)
    return T(mapping.domain, map(f -> refine_imp(f, method), mapping.data)...)
end

function refine_imp end
