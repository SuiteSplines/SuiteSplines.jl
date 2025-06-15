export Accessor, ElementAccessor, Element, Elements, TrialIndices, TestIndices, TrialFunctions, TestFunctions, QuadraturePoints, QuadratureWeights
export get_element_domain, get_element_interval

abstract type Accessor{Dim} end

"""
    ElementAccessor{Dim, X<:IncreasingSequence, Data}

Datastructure that provides access to element data
such as the element domain, element trialfunctions,
element testfunctions and element quadrature rule.
"""
struct ElementAccessor{Dim,X,U,CInd,LInd,Data} <:Accessor{Dim}
    partition::X        # partition
    testfuns::U         # test function space
    trialfuns::U        # trial function space
    eind::CInd          # element indices
    uind::LInd          # trialfunction indices
    vind::LInd          # testfunction indices
    data::Data          # Univariate data
    function ElementAccessor(x::X, v::U, u::U, ei::CInd, ui::LInd, vi::LInd, data::Data) where {X,U,LInd,CInd,Data}
        Dim = ndims(x)
        @assert ndims(u) == ndims(v) == ndims(data) == Dim
        @assert ndims(ui) == ndims(vi) == ndims(ei) == Dim
        return new{Dim,X,U,CInd,LInd,Data}(x, v, u, ei, ui, vi, data)
    end
end

Base.ndims(acc::ElementAccessor) = ndims(acc.partition)
Base.length(acc::ElementAccessor) = length(acc.eind)
Base.size(acc::ElementAccessor) = size(acc.eind)
Base.size(acc::ElementAccessor, k) = size(acc.eind, k)

function ElementAccessor(; testspace, trialspace, quadrule, kwargs...)
    partition = CartesianProduct(v -> breakpoints(v), trialspace) # ToDo check equality with breakpoints of testspace
    Dim = ndims(partition)
    eind = CartesianIndices(ntuple(k -> Base.OneTo(size(partition,k)-1), Dim))
    uind = LinearIndices(ntuple(i -> Base.OneTo(size(trialspace[i].C,1)), Dim))
    vind = LinearIndices(ntuple(i -> Base.OneTo(size(testspace[i].C,1)), Dim))
    data = ElementAccessorData(partition, trialspace, testspace, quadrule; kwargs...) # dispatch for different types of quadrature rules
    return ElementAccessor(partition, testspace, trialspace, eind, uind, vind, data)
end

"""
    Element{Dim, T}

Dim-dimensional element storing its domain of definition as a
`Cartesian{Dim, Interval{T}}` and its element number as a 
`CartesianIndex{Dim}`
"""
struct Element{Dim, X<:AbstractArray}
    index::CartesianIndex{Dim}
    parent::X
    function Element(index::CartesianIndex{Dim}, parent::AbstractArray{T,Dim}) where {T,Dim}
        X = typeof(parent)
        return new{Dim,X}(index, parent)
    end
end

function Element(acc::ElementAccessor, lindex::Int)
    indices = CartesianIndices(acc.partition)[lindex]
    return Element(indices, acc.partition)
end

Base.ndims(e::Element) = Base.ndims(e.parent)

@eval @inline function Base.getindex(element::Element{1}, i::Int)
    @assert i in 1:2
    return element.parent[element.index[1]-1+i]
end

@eval @inline function Base.getindex(element::Element{2}, i::Vararg{Int, 2})
    @assert i[1] in 1:2 && i[2] in 1:2
    return element.parent[element.index[1]+i[1]-1,element.index[2]+i[2]-1]
end

@eval @inline function Base.getindex(element::Element{3}, i::Vararg{Int, 3})
    @assert i[1] in 1:2 && i[2] in 1:2 && i[3] in 1:2
    return element.parent[element.index[1]+i[1]-1,element.index[2]+i[2]-1, element.index[3]+i[3]-1]
end

function Base.show(io::IO, element::Element)
    I = element.index
    parent = element.parent
    indices = ntuple(k -> I[k]:I[k]+1, ndims(element))
    corner_pts = parent[indices...]
    print(io, "Element with corner points: \n")
    Base.print_array(io, corner_pts)
end

function get_element_interval(element::Element{1}, dir)
    a = element[1]
    b = element[2]
    return Interval(a[dir], b[dir])
end

function get_element_interval(element::Element{2}, dir)
    a = element[1,1]
    b = element[2,2]
    return Interval(a[dir], b[dir])
end

function get_element_interval(element::Element{3}, dir)
    a = element[1,1,1]
    b = element[2,2,2]
    return Interval(a[dir], b[dir])
end

function get_element_domain(element)
    return CartesianProduct(dir -> get_element_interval(element,dir), 1:ndims(element))
end

"""
    coface(element::Element)

If `element` is on the boundary then `element.parent` is a view and
this function returns its single co-face.
"""
function coface(element::Element)
    parent = element.parent.parent
    indices = Base.reindex(element.parent.indices, element.index.I)
    indices = map((i,s) -> i - (i==s), indices, size(parent))
    return Element(CartesianIndex(indices), parent)
end

"""
    TrialIndices(acc::ElementAccessor{Dim}, element::Element{Dim})

Get global indices of the trial functions on an element.
"""
function TrialIndices(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim
    indices = ntuple(k -> get_basis_indices(acc.trialfuns[k], element.index[k]), Dim)
    return reshape(view(acc.uind, indices...), :)
end

TrialIndices(acc, element) = TrialIndices(acc, coface(element))

"""
    TestIndices(acc::ElementAccessor{Dim}, element::Element{Dim})

Get global indices of the test functions on an element.
"""
function TestIndices(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim
    indices = ntuple(k -> get_basis_indices(acc.testfuns[k], element.index[k]), Dim)
    return reshape(view(acc.vind, indices...), :)
end

TestIndices(acc, element) = TestIndices(acc, coface(element))

"""
    TrialFunctions(acc::ElementAccessor{Dim}, element::Element{Dim})

Get the trial functions on an element as a KroneckerProduct matrix.
"""
@inline function TrialFunctions(acc::ElementAccessor, e::Element; ders)
    return TrialFunctions(acc.data, e; ders=ders)
end

@inline function TrialFunctions(acc::ElementAccessor, e::Element, i::Int; ders)
    return TrialFunctions(acc.data, e, i; ders=ders)
end

"""
    TestFunctions(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim

Get the test functions on an element as a KroneckerProduct matrix.
"""
@inline function TestFunctions(acc::ElementAccessor, e::Element; ders)
    return TestFunctions(acc.data, e; ders=ders)
end

@inline function TestFunctions(acc::ElementAccessor, e::Element, i::Int; ders)
    return TestFunctions(acc.data, e, i; ders=ders)
end

"""
    QuadratureRule(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim

Get the quadrature rule on this element.
"""
@inline function IgaBase.QuadratureRule(acc::ElementAccessor, e::Element; kwargs...)
    return QuadratureRule(acc.data, e; kwargs...)
end

"""
    QuadraturePoints(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim

Get the quadrature points on this element.
"""
@inline function QuadraturePoints(acc::ElementAccessor, e::Element; kwargs...)
    return QuadraturePoints(acc.data, e; kwargs...)
end

"""
    QuadratureWeights(acc::ElementAccessor{Dim}, element::Element{Dim}) where Dim

Get the quadrature weights on this element.
"""
@inline function QuadratureWeights(acc::ElementAccessor, e::Element; kwargs...)
    return QuadratureWeights(acc.data, e; kwargs...)
end

"""
    Elements(partition::CartesianProduct{Dim, <:IncreasingSequence{<:Real}})

Return an iterator over all elements in a partition.
"""
struct Elements{Dim, X, Indices}
    partition::X
    eind::Indices
    function Elements(partition::X, eind::Indices) where {X, Indices}
        @assert ndims(partition) == ndims(eind)
        Dim = ndims(partition)
        new{Dim,X,Indices}(partition, eind)
    end
end

function Elements(partition::AbstractArray)
    return Elements(partition, CartesianIndices(ntuple(k -> Base.OneTo(max(1,size(partition, k)-1)), ndims(partition))))
end

Base.length(E::Elements) = length(E.eind)
Base.size(E::Elements) = size(E.eind)
Base.size(E::Elements, k) = size(E.eind, k)
Base.eltype(E::Elements{Dim}) where Dim = Element{Dim, eltype(E.partition[1])}

Base.checkbounds(E::Elements, I) = Base.checkbounds(E.eind, I...)

function Base.getindex(E::Elements, k::Int)
    @boundscheck checkbounds(E, k)
    return Element(E.eind[k], E.partition)
end

function Base.getindex(E::Elements{Dim}, I::Vararg{Int,Dim}) where Dim
    @boundscheck checkbounds(E, I)
    return Element(E.eind[I...], E.partition)
end

function Base.iterate(iter::Elements)
    return (iter[1], 1)
end

function Base.iterate(iter::Elements, state)
    while state<length(iter)
        state += 1
        return (iter[state], state)
    end
end

function get_basis_indices(S::SplineSpace, e::Integer)
    span = S.s[e]
    return span-Degree(S):span
end