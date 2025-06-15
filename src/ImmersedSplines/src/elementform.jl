export InsideOutsideTest, is_inside, is_outside, on_interface, CutcellQuadratureRule

struct InsideOutsideTest
    value::Int8
end

is_inside(test::InsideOutsideTest) = test.value==1
is_outside(test::InsideOutsideTest) = test.value==-1
on_interface(test::InsideOutsideTest) = test.value==0

# ToDo: Need to test at more points
function InsideOutsideTest(phi, element::Element{1})
    d = phi(SVector(element[1]))
    if d<0
        if phi(SVector(element[2])) > 0
            return InsideOutsideTest(0) # case on interface
        end
        return InsideOutsideTest(1) # case inside
    else
        if phi(SVector(element[2])) < 0
            return InsideOutsideTest(0) # case on interface
        end
        return InsideOutsideTest(-1) # case outside
    end
end

# ToDo: Need to test at more points
function InsideOutsideTest(phi, element::Element{2})
    d = phi(SVector(element[1,1]))
    if d<0
        for j in 1:2
            for i in 1:2
                if phi(SVector(element[i,j])) > 0
                    return InsideOutsideTest(0) # case on interface
                end
            end
        end
        return InsideOutsideTest(1) # case inside
    else
        for j in 1:2
            for i in 1:2
                if phi(SVector(element[i,j])) < 0
                    return InsideOutsideTest(0) # case on interface
                end
            end
        end
        return InsideOutsideTest(-1) # case outside
    end
end

# ToDo: Need to test at more points
function InsideOutsideTest(phi, element::Element{3})
    d = phi(SVector(element[1,1,1]))
    if d<0
        for k in 1:2
            for j in 1:2
                for i in 1:2
                    if phi(SVector(element[i,j,k])) > 0
                        return InsideOutsideTest(0) # case on interface
                    end
                end
            end
        end
        return InsideOutsideTest(1) # case inside
    else
        for k in 1:2
            for j in 1:2
                for i in 1:2
                    if phi(SVector(element[i,j,k])) < 0
                        return InsideOutsideTest(0) # case on interface
                    end
                end
            end
        end
        return InsideOutsideTest(-1) # case outside
    end
end

mutable struct CutcellQuadratureRule{Dim}
    partition::CartesianProduct{Dim}
    mapping::AlgoimCallLevelSetFunction
    npoints::Int
    element::Element
    quadrule::ImmersedQuadRule
    function CutcellQuadratureRule(; partition::CartesianProduct{Dim}, mapping, npoints) where {Dim}
        new{Dim}(partition, mapping, npoints)
    end
end

Base.isassigned(Q::CutcellQuadratureRule) = isdefined(Q, :element) && isdefined(Q, :quadrule)

function istodate(Q::CutcellQuadratureRule, e::Element, phase)::Bool
    
    # check if struct is already initialized
    if !isassigned(Q)
        return false
    end

    # check if the same quadrature rule is stored
    if Q.element==e && Q.quadrule.phase==phase
        return true
    end

    return false
end

function update!(Q::CutcellQuadratureRule, e::Element, phase::Int=-1)
    if !istodate(Q, e, phase)
        Q.element = e # element and quadrule are always updated simultaneously
        Q.quadrule = ImmersedQuadRule(Q.mapping, Q.element; order=Q.npoints, phase=phase)
    end
end

struct CutcellData{Dim, T, Q<:CutcellQuadratureRule}
    trialspace::TensorProduct{Dim,SplineSpace{T}}
    testspace::TensorProduct{Dim,SplineSpace{T}}
    quadrule::Q
end

Base.ndims(::CutcellData{Dim}) where Dim = Dim

function IgaFormation.ElementAccessorData(partition, trialspace, testspace, quadrule::CutcellQuadratureRule)
    return CutcellData(trialspace, testspace, quadrule)
end

function IgaFormation.QuadratureRule(data::CutcellData,  e::Element; phase::Int=-1)
    update!(data.quadrule, e, phase) # update algoim data 
    return data.quadrule.quadrule
end

function IgaFormation.QuadraturePoints(data::CutcellData,  e::Element; phase::Int=-1)
    Q = QuadratureRule(data, e; phase=phase)
    return Q.x
end

function IgaFormation.QuadratureWeights(data::CutcellData,  e::Element; phase::Int=-1)
    Q = QuadratureRule(data, e; phase=phase)
    return Q.w
end

function IgaFormation.TrialFunctions(data::CutcellData{Dim}, e::Element, i::Int; ders) where Dim
    x = single_quadrature_point(data, e, i) 
    eindex = IgaFormation.get_element_indices(e)
    funs = ntuple(k -> get_element_functions(data.trialspace[k], eindex[k], x[k], ders[k]), Dim)
    return KroneckerProduct(funs..., reverse=true)
end

function IgaFormation.TestFunctions(data::CutcellData{Dim}, e::Element, i::Int; ders) where Dim
    x = single_quadrature_point(data, e, i) 
    eindex = IgaFormation.get_element_indices(e)
    funs = ntuple(k -> get_element_functions(data.testspace[k], eindex[k], x[k], ders[k]), Dim)
    return KroneckerProduct(funs..., reverse=true)
end

function single_quadrature_point(data::CutcellData{Dim},  e::Element{Dim}, i::Int) where Dim
    @assert data.quadrule.element==e # check if data is up to date
    return data.quadrule.quadrule.x[i]
end

function single_quadrature_point(data::CutcellData{2},  e::Element{1}, i::Int)
    @assert data.quadrule.element==e # check if data is up to date
    k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
    x = data.quadrule.quadrule.x[i]
    s = e[1][k]
    if k==1
        return (s, x[1])
    elseif k==2
        return (x[1], s)
    end
end

function single_quadrature_point(data::CutcellData{3},  e::Element{2}, i::Int)
    @assert data.quadrule.element==e # check if data is up to date
    k = IgaBase.find_singular_dimension(map(length, e.parent.indices)...)
    x = data.quadrule.quadrule.x[i]
    s = e[1,1][k]
    if k==1
        return (s, x[1], x[2])
    elseif k==2
        return (x[1], s, x[2])
    elseif k==3
        return (x[1], x[2], s)
    end
end

function get_element_functions(S::SplineSpace, e::Int, u, ders::Int)
    return dersbsplinebasisfuns(Degree(S), KnotVector(S), u, ders+1)[:,ders+1]
end

# function get_span_indices(S::SplineSpace, e::Integer)
#     if e==0         # left boundary quadrature point
#         return S.s[1]
#     elseif e==length(S.s)
#         return S.s[end]
#     else    # interior quadrature points
#         @show e
#         @show length(S.s)
#         return S.s[e]
#     end
# end
