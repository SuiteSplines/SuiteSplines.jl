export PatchAccessor, TrialIndices, TestIndices, TrialFunctions, TestFunctions, QuadraturePoints, QuadratureWeights

"""
    PatchAccessor{Dim, X, ..., Data}

Datastructure that provides access to element data
such as the element domain, element trialfunctions,
element testfunctions and element quadrature rule.
"""
struct PatchAccessor{Dim,X,U,LInd,Data} <: Accessor{Dim}
    partition::X        # partition
    testfuns::U         # test function space
    trialfuns::U        # trial function space
    uind::LInd          # trialfunction indices
    vind::LInd          # testfunction indices
    data::Data          # Univariate data
    function PatchAccessor(x::X, v::U, u::U, ui::LInd, vi::LInd, data::Data) where {X,U,LInd,Data}
        Dim = ndims(x)
        @assert ndims(u) == ndims(v) == ndims(data) == Dim
        @assert ndims(ui) == ndims(vi)== Dim
        return new{Dim,X,U,LInd,Data}(x, v, u, ui, vi, data)
    end
end

Base.ndims(acc::PatchAccessor) = ndims(acc.partition)
Base.length(acc::PatchAccessor) = length(acc.eind)
Base.size(acc::PatchAccessor) = size(acc.eind)
Base.size(acc::PatchAccessor, k) = size(acc.eind, k)

function PatchAccessor(; testspace, trialspace, quadrule, kwargs...)
    partition = CartesianProduct(v -> breakpoints(v), trialspace) # ToDo check equality with breakpoints of testspace
    Dim = ndims(partition)
    uind = LinearIndices(ntuple(i -> Base.OneTo(dimsplinespace(trialspace[i])), Dim))
    vind = LinearIndices(ntuple(i -> Base.OneTo(dimsplinespace(testspace[i])), Dim))
    data = PatchAccessorData(partition, trialspace, testspace, quadrule; kwargs...) # dispatch for different types of quadrature rules
    return PatchAccessor(partition, testspace, trialspace, uind, vind, data)
end

"""
    TrialIndices(acc::PatchAccessor)

Get global indices of the trial functions.
"""
function TrialIndices(acc::ElementAccessor)
    return acc.uind
end

"""
    TestIndices(acc::PatchtAccessor)

Get global indices of the test functions.
"""
function TestIndices(acc::ElementAccessor)
    return acc.vind;
end


"""
    TrialFunctions(acc::PatchAccessor)

Get the trial functions on the patch as a KroneckerProduct matrix.
"""
@inline function TrialFunctions(acc::PatchAccessor; ders)
    return TrialFunctions(acc.data; ders=ders)
end

"""
    TestFunctions(acc::PatchAccessor{Dim}) where Dim

Get the test functions on the patch as a KroneckerProduct matrix.
"""
@inline function TestFunctions(acc::PatchAccessor; ders)
    return TestFunctions(acc.data; ders=ders)
end

"""
    QuadratureRule(acc::PatchAccessor{Dim}) where Dim

Get the quadrature rule on the patch.
"""
@inline function IgaBase.QuadratureRule(acc::PatchAccessor)
    return QuadratureRule(acc.data)
end

"""
    QuadraturePoints(acc::PatchAccessor{Dim}) where Dim

Get the quadrature points on this element.
"""
@inline function QuadraturePoints(acc::PatchAccessor)
    return QuadraturePoints(acc.data)
end

"""
    QuadratureWeights(acc::PatchAccessor{Dim}) where Dim

Get the quadrature weights on this element.
"""
@inline function QuadratureWeights(acc::PatchAccessor)
    return QuadratureWeights(acc.data)
end