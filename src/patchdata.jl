export TrialFunctions, TestFunctions, QuadraturePoints, QuadratureWeights, QuadratureRule

function PatchAccessorData(partition, trialspace, testspace, quadrule::TensorProduct; kwargs...)
    return TensorProduct((u, Su, Sv, qr) -> UnivariatePatchData(u, Su, Sv, qr; kwargs...), partition, trialspace, testspace, quadrule)
end

function TrialFunctions(data::CartesianPatchData{Dim}; ders::NTuple{Dim,Int}) where Dim
    return KroneckerProduct(k -> get_trialbasis(data[k], ders[k]), 1:Dim; reverse=true)
end

function TestFunctions(data::CartesianPatchData{Dim}; ders::NTuple{Dim,Int}) where Dim
    return KroneckerProduct(k -> get_testbasis(data[k], ders[k]), 1:Dim; reverse=true)
end

function QuadraturePoints(data::CartesianPatchData)
    return CartesianProduct(k -> get_qpoints(data[k]), 1:ndims(data))
end

function QuadratureWeights(data::CartesianPatchData)
    return KroneckerProduct(k -> get_qweights(data[k]), 1:ndims(data); reverse=true)
end

function IgaBase.QuadratureRule(data::CartesianPatchData)
    return QuadratureRule(QuadraturePoints(data), QuadratureWeights(data))
end

function get_testspace_dimension(patchdata::UnivariatePatchData)
    return dimsplinespace(patchdata.testspace)
end

function get_trialspace_dimension(patchdata::UnivariatePatchData)
    return dimsplinespace(patchdata.trialspace)
end

function get_n_quadpoints(patchdata::UnivariatePatchData)
    return length(patchdata.qrule.x)-2
end

function get_trialbasis(patchdata::UnivariatePatchData, ders::Int=0)
    return Matrix(patchdata.trialspace.C' * patchdata.trialbasis[:,2:end-1,ders+1])
end

function get_testbasis(patchdata::UnivariatePatchData, ders::Int=0)
    return Matrix(patchdata.testspace.C' * patchdata.testbasis[:,2:end-1,ders+1])
end

function get_qpoints(patchdata::UnivariatePatchData)
    return @view patchdata.qrule.x[2:end-1] 
end

function get_qweights(patchdata::UnivariatePatchData)
    return @view patchdata.qrule.w[2:end-1]
end