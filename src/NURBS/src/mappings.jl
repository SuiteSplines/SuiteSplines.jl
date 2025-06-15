export Nurbs, transform_to_projective_space!, Mapping, Boundary

# NURBS:
# This does nothing for NURBS.. why is it here?
#function AbstractMappings.check_mapping_arguments(b::NTuple{Codim,TensorProductBspline}) where Codim
#    cache = b[1].cache
#    space = b[1].space
#    for k in 1:Codim
#        @assert sum(b[k].ders) == 0 "Provide B-spline functions; not their derivatives."
#        @assert b[k].cache === cache "Provide B-splines with a single cache."
#        @assert b[k].space === space "Provide B-splines defined via the same space."
#    end
#end

function AbstractMappings.get_property_imp(::Type{<:NTuple{Codim,Nurbs}}, f::GeometricMapping, s::Symbol) where Codim
    s === :space        && return getproperty(getfield(f, :data)[1], :space)
    s === :weightfun    && return getfield(getfield(f, :data)[1], :weightfun)
    s === :splinefun    && return getfield(getfield(f, :data)[1], :splinefun)
    s === :cache        && return getproperty(getfield(f, :data)[1], :cache)
    s === :weights      && return getproperty(getfield(f, :data)[1], :weights)
    getfield(f, s)
end

function transform_to_projective_space!(F::AbstractMapping)
    for f in F
        transform_to_projective_space!(f)
    end
end

# Boundary iterator
function IgaBase.boundary_element_type(s::Nurbs{Dim}) where Dim
    S = IgaBase.boundary_element_type(s.splinefun)
    return Nurbs{Dim-1,S}
end