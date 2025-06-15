export ElementAccessor, Element, Elements, TrialIndices, TestIndices, TrialFunctions, TestFunctions, QuadraturePoints, QuadratureWeights

function ElementAccessorData(partition, trialspace, testspace, quadrule::TensorProduct; kwargs...)
    return TensorProduct((u, Su, Sv, qr) -> UnivariatePatchData(u, Su, Sv, qr; kwargs...), partition, trialspace, testspace, quadrule)
end

"""
    UnivariatePatchData{T<:Real, Q<:AbstractQuadrule{1}}

Datastructure that stores all the univariate patch data. This data
is conveniently accessed by an element accessor.
"""
struct UnivariatePatchData{T<:Real, X<:IncreasingSequence{T}, Q<:AbstractQuadrule{1}}
    partition::X
    testspace::SplineSpace{T}
    trialspace::SplineSpace{T}
    testbasis::BsplineBasis{T}
    trialbasis::BsplineBasis{T}
    qrule::Q
end

Base.ndims(::UnivariatePatchData) = 1

function UnivariatePatchData(partition, trialspace, testspace, quadrule::AbstractQuadrule{1}; incorporate_weights_in_testfuns=false)
    trialbasis = BsplineBasis(trialspace, quadrule.x, trialspace.p+1)
    testbasis  = BsplineBasis(testspace, quadrule.x, testspace.p+1)
    if (incorporate_weights_in_testfuns)
        testbasis.data[:,2:end-1,:] .*= quadrule.w[2:end-1]'
    end
    UnivariatePatchData(partition, testspace, trialspace, testbasis, trialbasis, quadrule)
end

const CartesianPatchData{Dim,T,X,Q} = TensorProduct{Dim,UnivariatePatchData{T,X,Q}}

function TrialFunctions(data::CartesianPatchData{Dim}, element::Element; ders::NTuple{Dim,Int}) where Dim
    eindex = get_element_indices(element)
    return KroneckerProduct(k -> get_element_trialbasis(data[k], eindex[k], ders[k]), 1:Dim; reverse=true)
end

function TestFunctions(data::CartesianPatchData{Dim}, element::Element; ders::NTuple{Dim,Int}) where Dim
    eindex = get_element_indices(element)
    return KroneckerProduct(k -> get_element_testbasis(data[k], eindex[k], ders[k]), 1:Dim; reverse=true)
end

function QuadraturePoints(data::CartesianPatchData, element::Element; kwargs...)
    eindex = get_element_indices(element)
    return CartesianProduct(k -> get_element_qpoints(data[k], eindex[k]), 1:ndims(data))
end

function QuadratureWeights(data::CartesianPatchData, element::Element; kwargs...)
    eindex = get_element_indices(element)
    return KroneckerProduct(k -> get_element_qweights(data[k], eindex[k]), 1:ndims(data); reverse=true)
end

function IgaBase.QuadratureRule(data::CartesianPatchData, element::Element; kwargs...)
    return QuadratureRule(QuadraturePoints(data, element), QuadratureWeights(data, element))
end

function get_element_indices(element)
    return element.index.I
end

function get_element_indices(element::Element{1,<:SubArray})
    k = IgaBase.find_singular_dimension(map(length,element.parent.indices)...)
    s = element.parent.indices[k]; s-=(s==1)
    if k==1
        return (s, element.index[1])
    elseif k==2
        return (element.index[1], s)
    else
        return element.index.i
    end
end

function get_element_indices(element::Element{2,<:SubArray})
    k = IgaBase.find_singular_dimension(map(length,element.parent.indices)...)
    s = element.parent.indices[k]; s-=(s==1)
    if k==1
        return (s, element.index[1], element.index[2])
    elseif k==2
        return (element.index[1], s, element.index[2])
    elseif k==3
        return (element.index[1], element.index[2], s)
    else
        return element.index.i
    end
end

function get_n_elements(patchdata::UnivariatePatchData)
    return length(patchdata.partition)-1
end

function get_element_testspace_dimension(patchdata::UnivariatePatchData)
    return Degree(patchdata.testspace)+1
end

function get_element_trialspace_dimension(patchdata::UnivariatePatchData)
    return Degree(patchdata.trialspace)+1
end

function get_element_n_quadpoints(patchdata::UnivariatePatchData, e::Integer)
    @assert 0 <= e <= get_n_elements(patchdata)+1
    return patchdata.qrule.i[e+2] - patchdata.qrule.i[e+1]
end

function get_max_n_quadpoints(patchdata::UnivariatePatchData)
    nquad = 0;
    for e in 1:get_n_elements(patchdata)
        l = get_element_n_quadpoints(patchdata, e)
        if l>nquad
            nquad=l
        end
    end
    return nquad
end

function get_span_indices(patchdata::UnivariatePatchData, e::Integer)
    if e==0         # left boundary quadrature point
        return 1:1
    elseif e==get_n_elements(patchdata)+1   # right boundary quadrature point
        return length(patchdata.qrule.x):length(patchdata.qrule.x)
    else    # interior quadrature points
        return patchdata.qrule.i[e+1]:patchdata.qrule.i[e+2]-1 # get span indices of quadrature points
    end
end

function get_testfunction_spanindex(patchdata::UnivariatePatchData, e::Integer)
    @assert 0 < e <= get_n_elements(patchdata) # check bounds of element number
    p = patchdata.testspace.p
    span = patchdata.testspace.s[e]
    return span-p:span
end

function get_trialfunction_spanindex(patchdata::UnivariatePatchData, e::Integer)
    @assert 0 < e <= get_n_elements(patchdata) # check bounds of element number
    p = patchdata.trialspace.p
    span = patchdata.trialspace.s[e]
    return span-p:span
end

function get_element_interval(patchdata::UnivariatePatchData, e::Integer)
    @assert 0 < e <= get_n_elements(patchdata) # check bounds of element number
    return Interval(patchdata.partition[e], patchdata.partition[e+1]) # get element interval
end

function get_element_trialbasis(patchdata::UnivariatePatchData, e::Integer, ders::Int=0)
    indices = get_span_indices(patchdata, e)
    return view(patchdata.trialbasis.data, :, indices, ders+1) # element trial functions
end

function get_element_testbasis(patchdata::UnivariatePatchData, e::Integer, ders::Int=0)
    indices = get_span_indices(patchdata, e)
    return view(patchdata.testbasis.data, :, indices, ders+1) # element trial functions
end

function get_element_qpoints(patchdata::UnivariatePatchData, e::Integer)
    indices = get_span_indices(patchdata, e)
    return view(patchdata.qrule.x, indices) # element quadrature points 
end

function get_element_qweights(patchdata::UnivariatePatchData, e::Integer)
    if e==0 || e==get_n_elements(patchdata)+1
        return [1.0]
    else
        return view(patchdata.qrule.w, get_span_indices(patchdata, e)) # element quadrature points 
    end
end