
Base.convert(::Type{S}, v::CartesianProduct{1,T,Tuple{S}}) where {T,S<:IncreasingSequence{T}} = v.data[1] 
Base.convert(::Type{T}, v::TensorProduct{1,T}) where {T} = v.data[1]  
