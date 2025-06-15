export Pattern, kformbasisconv

using Base.Cartesian

"""
    subsets(v::NTuple{n,Int}, k::Int)

Generate the `binomial(n,k)` unique subsets of length `k` from set `v` and its
complement.

# Examples:
```@jldoctest
julia> α, β = subsets((1,2,3), 1);

julia> α
3-element Array{Tuple{Int64},1}:
 (1,)
 (2,)
 (3,)

julia> β
3-element Array{Tuple{Int64,Int64},1}:
 (2, 3)
 (1, 3)
 (1, 2)
```
"""
@generated function subsets(v::NTuple{N,Int}, ::Val{K}) where {N,K}
    quote
        @nexprs 1 j->(i_{j+$K}=0)
        Σ = NTuple{$K,Int}[]
        Τ = NTuple{$N-$K,Int}[]
        @nloops $K i j->i_{j+1}+1:$N+1-j begin
            tmp = @ntuple $K i -> v[i_{$K+1-i}]
            push!(Σ, tmp)
            push!(Τ, tuple(setdiff(v,tmp)...))
        end
        return Σ, Τ
    end
end

subsets(v::NTuple{N,Int}, k::Int) where N = subsets(v, Val(k))

supp(a) = [0:length(a)-1...][a.!=0]

# Defintion of Pattern - a short linear combination
struct Pattern{T<:Real}
    m::Integer
    B::Vector{Int}
    S::Vector{T}
    G::Array{Int}
end

"""
    kformbasisconv(n::Dimension, k::Form, r::Degree)

Compute the linear mapping from Bernstein polynomials to the space of
`k`-forms `Pᵣ⁻Λᵏ(T)`. This linear mapping is encoded in `Pattern{T}`. We
refer to the [paper](https://cpb-us-w2.wpmucdn.com/sites.baylor.edu/dist/6/32/files/2015/01/ssec-1zn3ii2.pdf).
For more information, check out the [periodic table of finite elements](http://femtable.org/background.html)
"""
function kformbasisconv(n::Dimension, k::Form, r::Degree)
    # compute the conversion matrix of the differential k-form to Bernstein polynomials
    # of degree r on the parent triangle
    # n - spatial dim
    # k - degree differential form
    # r - polynomial degree

    # determine the dim of the k-form space Pᵣ⁻Λᵏ(T)
    m = (k+1) * Binomial(r+k-1, k) * Binomial(n+r, n-k)

    # allocate space for vectors
    B = zeros(Int, m)
    S = zeros(Float64, m)
    G = zeros(Int, m, k)

    # canonical directions
    δ = collect(MultiIndex, MultiIndices(n+1))

    # determine face sets
    Σ, Τ = subsets(ntuple(i -> i-1, n+1), k+1)

    # determine pattern
    it = 1
    for σ in Iterators.map(σ -> [σ...], Σ)

        for α in MultiIndices(r-1,n+1)
            # test condition
            if all(α[1:minimum(σ)].==0)

                for i in 0:k
                    # counters
                    ii = i+1
                    σi = σ[ii]+1
                    j = (k+1)*(it-1) + i+1

                    # fill data structures
                    S[j] = (-1)^(i) * (α[σi] + 1) / r
                    B[j] = linindex(α .+ δ[σi])
                    G[j,:] .= σ[0:k .!= i]
                end

                it+=1
            end
        end
    end
    return Pattern(k+1,B,S,G)
end
