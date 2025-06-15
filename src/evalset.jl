export EvaluationSet, elsize

"""
EvaluationSet{S1,S2, S}

Set of evaluation points that are used as input and output of functions,
mappings and fields.
"""
struct EvaluationSet{S1,S2,T,Dims,Array<:AbstractArray} <: AbstractArray{T,Dims}
    data::SizedMatrix{S1,S2,Array}
    function EvaluationSet{S1,S2}(x::NTuple{L,Array}) where {S1,S2,L,Array}
        @assert L==S1*S2
        for k in 2:length(x)
            @assert size(x[k])==size(x[1])
        end 
        Dims = ndims(x[1])
        Tscalar = eltype(x[1])
        T = (S1==1 && S2==1) ? Tscalar : SMatrix{S1,S2,Tscalar}
        return new{S1,S2,T,Dims,Array}(SizedMatrix{S1,S2,Array}(x))
    end
end 

EvaluationSet{S1,S2}(x::Array) where {S1,S2,Array<:AbstractArray} =  EvaluationSet{S1,S2}((x,))
EvaluationSet{S1,S2}(x::Array...) where {S1,S2,Array<:AbstractArray} =  EvaluationSet{S1,S2}(x)

Base.eltype(Y::EvaluationSet{1,1}) = eltype(Y.data[1])
Base.eltype(Y::EvaluationSet{S1,S2}) where {S1,S2} = SMatrix{S1,S2,eltype(Y.data[1])}

import Base.elsize
Base.elsize(::EvaluationSet{S1,S2}) where {S1,S2} = (S1,S2)
Base.elsize(Y::EvaluationSet, k) = elsize(Y)[k]

Base.ndims(Y::EvaluationSet) = ndims(Y.data[1]) 
Base.length(Y::EvaluationSet) = length(Y.data[1])
Base.size(Y::EvaluationSet, k) = size(Y.data[1], k)
Base.size(Y::EvaluationSet) = size(Y.data[1])

EvaluationSet{S1,S2}(array_type, array_initializer, ndims...) where {S1,S2} = EvaluationSet{S1,S2}(ntuple(k -> array_type(array_initializer, ndims), S1*S2))

# getindex and setindex! methods
for S1 in 1:3
    for S2 in 1:3
        # special case scalar evaluation
        if S2==1 && S1==1
            @eval function Base.getindex(Y::EvaluationSet{1,1}, I...)
                return Y.data[1][I...]
            end
            @eval function Base.setindex!(Y::EvaluationSet{1,1}, v, I...)
                Y.data[1][I...] = v
            end
        # general case
        else        
            @eval function Base.getindex(Y::EvaluationSet{$S1,$S2}, I...)
                return @SMatrix [Y.data[i,j][I...] for i in 1:$S1, j in 1:$S2]
            end
            @eval function Base.setindex!(Y::EvaluationSet{$S1,$S2}, v, I...)
                for j in 1:$S2
                    for i in 1:$S1
                        Y.data[i,j][I...] = v[i,j]
                    end
                end
            end
        end
    end
end

extract_at(A::AbstractArray, i...) = getindex(A, i...)
extract_at(A::EvaluationSet, i...) = getindex(A.data, i...)