export sumfact!, extract_element_cache, extract_patch_cache, pre_allocate_arrays, pre_allocate_rhs_arrays, FieldEvaluationCache
export VectorSumfactoryCache, MatrixSumfactoryCache, VectorSumfactory, MatrixSumfactory, Sumfactory

struct FieldEvaluationCache{S<:Accessor, T<:AbstractArray}
    acc::S
    data::T
end

function FieldEvaluationCache(acc::ElementAccessor{Dim}, M::Int=1, N::Int=1) where Dim
    T = Float64
    q = ntuple(k -> IgaFormation.get_max_n_quadpoints(acc.data[k]), Dim)
    data = EvaluationSet{M,N}(ResizableArray{Float64,Dim}, undef, q...)
    return FieldEvaluationCache(acc, data)
end

function FieldEvaluationCache(acc::PatchAccessor{Dim}, M::Int=1, N::Int=1) where Dim
    q = ntuple(k -> IgaFormation.get_n_quadpoints(acc.data[k]), Dim)
    data = EvaluationSet{M,N}(ResizableArray{Float64,Dim}, undef, q...)
    return FieldEvaluationCache(acc, data)
end

function extract_element_cache(cache::FieldEvaluationCache, element::Element)
    eindex = get_element_indices(element) 
    q = ntuple(k -> get_element_n_quadpoints(cache.acc.data[k], eindex[k]), length(eindex))
    for k in eachindex(cache.data.data)
        resize!(cache.data.data[k], q)
    end
    return cache.data
end

function extract_patch_cache(cache::FieldEvaluationCache)
    return cache.data
end

struct ElementIterationTrait
end

struct PatchIterationTrait
end

abstract type SumfactoryCache{Dim} end

Base.ndims(cache::SumfactoryCache) = Base.ndims(cache.acc) 

struct MatrixSumfactoryCache{Dim, T<:Real, A<:Accessor{Dim}} <: SumfactoryCache{Dim}
    acc::A
    data::NTuple{Dim, Vector{T}}
    data_b::NTuple{Dim, Vector{T}}
end

struct VectorSumfactoryCache{Dim, T<:Real, A<:Accessor{Dim}} <: SumfactoryCache{Dim}
    acc::A
    data::NTuple{Dim, Vector{T}}
end

function MatrixSumfactoryCache(acc::ElementAccessor)

    T = Float64
    Dim = ndims(acc)

    q = ntuple(k -> IgaFormation.get_max_n_quadpoints(acc.data[k]), Dim)
    n = ntuple(k -> IgaFormation.get_element_trialspace_dimension(acc.data[k]), Dim)
    m = ntuple(k -> IgaFormation.get_element_testspace_dimension(acc.data[k]), Dim)

    a, b = get_sumfact_alloc_dims(m, n, q)

    A = ntuple(k -> zeros(T, prod(a[k])), Dim)
    B = ntuple(k -> zeros(T, prod(b[k])), Dim)
    
    return MatrixSumfactoryCache(acc, A, B)
end

function VectorSumfactoryCache(acc::ElementAccessor)

    T = Float64
    Dim = ndims(acc)

    q = ntuple(k -> IgaFormation.get_max_n_quadpoints(acc.data[k]), Dim)
    m = ntuple(k -> IgaFormation.get_element_testspace_dimension(acc.data[k]), Dim)

    f = get_sumfact_alloc_dims(m, q)
    
    F = ntuple(k -> zeros(T, prod(f[k])), Dim)
    
    return VectorSumfactoryCache(acc, F)
end

function MatrixSumfactoryCache(acc::PatchAccessor)

    T = Float64
    Dim = ndims(acc)

    q = ntuple(k -> IgaFormation.get_n_quadpoints(acc.data[k]), Dim)
    n = ntuple(k -> IgaFormation.get_trialspace_dimension(acc.data[k]), Dim)
    m = ntuple(k -> IgaFormation.get_testspace_dimension(acc.data[k]), Dim)

    a, b = get_sumfact_alloc_dims(m, n, q)

    A = ntuple(k -> zeros(T, prod(a[k])), Dim)
    B = ntuple(k -> zeros(T, prod(b[k])), Dim)
    
    return MatrixSumfactoryCache(acc, A, B)
end

function VectorSumfactoryCache(acc::PatchAccessor)

    T = Float64
    Dim = ndims(acc)

    q = ntuple(k -> IgaFormation.get_n_quadpoints(acc.data[k]), Dim)
    m = ntuple(k -> IgaFormation.get_testspace_dimension(acc.data[k]), Dim)

    f = get_sumfact_alloc_dims(m, q)
    F = ntuple(k -> zeros(T, prod(f[k])), Dim)
    
    return VectorSumfactoryCache(acc, F)
end

abstract type AbstractSumfactory{Dim} end

struct MatrixSumfactory{Dim, T, U<:NTuple{Dim,Array{T}}, V<:NTuple{Dim, Matrix{T}}} <: AbstractSumfactory{Dim}
    data::U
    data_b::V
end

struct VectorSumfactory{Dim, T, U<:NTuple{Dim,Array{T}}} <: AbstractSumfactory{Dim}
    data::U
end

Base.ndims(sumfact::AbstractSumfactory) = length(sumfact.data)

allows_element_iteration(::ElementAccessor) = ElementIterationTrait()
allows_element_iteration(cache::SumfactoryCache) = allows_element_iteration(cache.acc)

allows_patch_iteration(::PatchAccessor) = PatchIterationTrait()
allows_patch_iteration(cache::SumfactoryCache) = allows_patch_iteration(cache.acc)


function MatrixSumfactory(::ElementIterationTrait, cache::MatrixSumfactoryCache, element::Element)
    
    # initialize
    Dim = ndims(cache)

    # reset element matrices and vectors
    reset!(cache)

    # get element particulars
    eindex = IgaFormation.get_element_indices(element)
    q = ntuple(k -> get_element_n_quadpoints(cache.acc.data[k], eindex[k]), Dim)
    n = ntuple(k -> get_element_trialspace_dimension(cache.acc.data[k]), Dim)
    m = ntuple(k -> get_element_testspace_dimension(cache.acc.data[k]), Dim)

    a, b = get_sumfact_alloc_dims(m, n, q)

    data  = ntuple(k -> resize!(cache.data[k], a[k]), Dim)
    data_b  = ntuple(k -> resize!(cache.data_b[k], b[k]), Dim)
    
    return MatrixSumfactory(data, data_b)
end

function VectorSumfactory(::ElementIterationTrait, cache::VectorSumfactoryCache, element::Element)
    
    # initialize
    Dim = ndims(cache)

    # reset element matrices and vectors
    reset!(cache)

    # get element particulars
    eindex = IgaFormation.get_element_indices(element)
    q = ntuple(k -> get_element_n_quadpoints(cache.acc.data[k], eindex[k]), Dim)
    m = ntuple(k -> get_element_testspace_dimension(cache.acc.data[k]), Dim)

    a = get_sumfact_alloc_dims(m, q)

    data  = ntuple(k -> resize!(cache.data[k], a[k]), Dim)
    
    return VectorSumfactory(data)
end

Sumfactory(cache::MatrixSumfactoryCache, element::Element) = MatrixSumfactory(allows_element_iteration(cache), cache, element)
Sumfactory(cache::VectorSumfactoryCache, element::Element) = VectorSumfactory(allows_element_iteration(cache), cache, element)



function MatrixSumfactory(::PatchIterationTrait, cache::MatrixSumfactoryCache)
    
    # initialize
    Dim = ndims(cache)

    # reset element matrices and vectors
    reset!(cache)

    # get element particulars
    q = ntuple(k -> get_n_quadpoints(cache.acc.data[k]), Dim)
    n = ntuple(k -> get_trialspace_dimension(cache.acc.data[k]), Dim)
    m = ntuple(k -> get_testspace_dimension(cache.acc.data[k]), Dim)

    a, b = get_sumfact_alloc_dims(m, n, q)

    data  = ntuple(k -> resize!(cache.data[k], a[k]), Dim)
    data_b  = ntuple(k -> resize!(cache.data_b[k], b[k]), Dim)
    
    return MatrixSumfactory(data, data_b)
end

function VectorSumfactory(::PatchIterationTrait, cache::VectorSumfactoryCache)
    
    # initialize
    Dim = ndims(cache)

    # reset element matrices and vectors
    reset!(cache)

    # get element particulars
    q = ntuple(k -> get_n_quadpoints(cache.acc.data[k]), Dim)
    m = ntuple(k -> get_testspace_dimension(cache.acc.data[k]), Dim)

    a = get_sumfact_alloc_dims(m, q)

    data  = ntuple(k -> resize!(cache.data[k], a[k]), Dim)
    
    return VectorSumfactory(data)
end

Sumfactory(cache::MatrixSumfactoryCache) = MatrixSumfactory(allows_patch_iteration(cache), cache)
Sumfactory(cache::VectorSumfactoryCache) = VectorSumfactory(allows_patch_iteration(cache), cache)



function reset!(cache::Union{MatrixSumfactory,MatrixSumfactoryCache}, hard_reset=true)
    # reset all except the output array of sum factorization
    for i in ndims(cache):-1:2
        cache.data[i] .= 0
        cache.data_b[i] .= 0
    end
    cache.data_b[1] .= 0

    # in a hard reset the output array of sum-factorization is also reset
    # if hard_reset==false then the previous results may be appended
    if hard_reset
        cache.data[1] .= 0
    end
end

function reset!(cache::Union{VectorSumfactory,VectorSumfactoryCache}, hard_reset=true)
    # reset all except the output array of sum factorization
    for i in ndims(cache):-1:2
        cache.data[i] .= 0
    end

    # in a hard reset the output array of sum-factorization is also reset
    # if hard_reset==false then the previous results may be appended
    if hard_reset
        cache.data[1] .= 0
    end
end

function (sumfact::MatrixSumfactory)(trialfuns::KroneckerProduct, testfuns::KroneckerProduct; data::DenseArray, reset=true)
    reset!(sumfact, reset)
    sumfact!(sumfact.data, sumfact.data_b, data, reverse(testfuns.data), reverse(trialfuns.data))
end

function (sumfact::VectorSumfactory)(testfuns::KroneckerProduct; data::DenseArray, reset=true)
    reset!(sumfact, reset)
    sumfact!(sumfact.data, data, reverse(testfuns.data))
end

function get_sumfact_alloc_dims(m::NTuple{3, Int}, n::NTuple{3, Int}, q::NTuple{3, Int})
    a₃ = (m[3], n[3], q[1], q[2])             # array sumfactorization stage 3
    a₂ = (m[3]*m[2], n[3]*n[2], q[1])         # array sumfactorization stage 2
    a₁ = (m[3]*m[2]*m[1], n[3]*n[2]*n[1])     # array sumfactorization stage 1
    b  = ntuple(k -> (m[k], n[k]), 3)         # container for outer products
    return (a₁, a₂, a₃), b
end

function get_sumfact_alloc_dims(m::NTuple{3, Int}, q::NTuple{3, Int})
    a₃ = (m[3], q[1], q[2])         # array sumfactorization stage 3
    a₂ = (m[3]*m[2], q[1])          # array sumfactorization stage 2
    a₁ = (m[3]*m[2]*m[1],)          # array sumfactorization stage 1
    return (a₁, a₂, a₃)
end

function get_sumfact_alloc_dims(m::NTuple{2, Int}, n::NTuple{2, Int}, q::NTuple{2, Int})
    a₂ = (m[2], n[2], q[1])             # array sumfactorization stage 2
    a₁ = (m[2]*m[1], n[2]*n[1])         # array sumfactorization stage 1
    b  = ntuple(k -> (m[k], n[k]), 2)   # container for outer products
    return (a₁, a₂), b
end

function get_sumfact_alloc_dims(m::NTuple{2, Int}, q::NTuple{2, Int})
    a₂ = (m[2], q[1])         # array sumfactorization stage 2
    a₁ = (m[2]*m[1],)          # array sumfactorization stage 1
    return (a₁, a₂)
end


"""
    pre_allocate_arrays(T, nquad::NTuple{3,Int}, m::NTuple{3,Int}, n::NTuple{3,Int})

Preallocate arrays of element type==T for use in sum factorization of matrices.
"""
function pre_allocate_arrays(T, nquad::NTuple{3,Int}, m::NTuple{3,Int}, n::NTuple{3,Int})
    A₃ = zeros(T, m[3], n[3], nquad[1], nquad[2])         # array sumfactorization stage 3
    A₂ = zeros(T, m[3]*m[2], n[3]*n[2], nquad[1])         # array sumfactorization stage 2
    A₁ = zeros(T, m[3]*m[2]*m[1], n[3]*n[2]*n[1])         # array sumfactorization stage 1
    B  = ntuple(k -> zeros(T, m[k], n[k]), 3)             # container for outer products
    return (A₁, A₂, A₃), B
end


"""
    pre_allocate_arrays(T, nquad::NTuple{3,Int}, m::NTuple{3,Int})

Preallocate arrays of element type==T for use in sum factorization of vectors.
"""
function pre_allocate_arrays(T, nquad::NTuple{3,Int}, m::NTuple{3,Int})
    A₃ = zeros(T, m[3], nquad[1], nquad[2])     # array sumfactorization stage 3
    A₂ = zeros(T, m[3]*m[2], nquad[1])          # array sumfactorization stage 2
    A₁ = zeros(T, m[3]*m[2]*m[1])               # array sumfactorization stage 1
    return (A₁, A₂, A₃)
end

"""
    pre_allocate_arrays(T, nquad::NTuple{2,Int}, m::NTuple{2,Int}, n::NTuple{2,Int})

Preallocate arrays of element type==T for use in sum factorization of matrices.
"""
function pre_allocate_arrays(T, nquad::NTuple{2,Int}, m::NTuple{2,Int}, n::NTuple{2,Int})
    A₂ = zeros(T, m[2], n[2], nquad[1])         # array sumfactorization stage 2
    A₁ = zeros(T, m[2]*m[1], n[2]*n[1])         # array sumfactorization stage 1
    B  = ntuple(k -> zeros(T, m[k], n[k]), 2)   # container for outer products
    return (A₁, A₂), B
end


"""
    pre_allocate_arrays(T, nquad::NTuple{2,Int}, m::NTuple{2,Int})

Preallocate arrays of element type==T for use in sum factorization of vectors.
"""
function pre_allocate_arrays(T, nquad::NTuple{2,Int}, m::NTuple{2,Int})
    A₂ = zeros(T, m[2], nquad[1])               # array sumfactorization stage 2
    A₁ = zeros(T, m[2]*m[1])                    # array sumfactorization stage 1
    return (A₁, A₂)
end