export Field, Pairing, ∘

struct Field{Dim,S1,S2,T} <: AbstractMapping{Dim,S1,S2}
    data::SizedMatrix{S1,S2,T}
    function Field{S1,S2}(args...) where {S1,S2}
        @assert length(args)==S1*S2 "Field dimensions inconsistent."
        y = map(process_mapping_input, args)
        Dim = n_input_args(y...)
        A = SizedMatrix{S1,S2}(y...)
        T = eltype(A)
        return new{Dim,S1,S2,T}(A)
    end
end

Field(args...) = Field{length(args),1}(args...)

function Field(T::Type, args::Tuple; codimension::Tuple{Int,Int}=(1,1))
    S1, S2 = codimension
    return Field{S1,S2}(map(T, args)...)
end

function Field(T::Type, args...; codimension::Tuple{Int,Int}=(1,1))
    x = T(args...)
    y = ntuple(k -> similar(x), prod(codimension)-1)
    S1, S2 = codimension
    return Field{S1,S2}(x, y...)
end

function Field(T::Type, args::Tuple{S1}; codimension::Tuple{Int,Int}=(1,1)) where S1
    @assert codimension[1]==S1
    S2 = codimension[2]
    return Field{S1,S2}(map(T, args)...)
end

struct Pairing{Dim1,Dim2,S1,S2,U<:GeometricMapping{Dim1,Dim2}, V<:AbstractMapping{Dim2,S1,S2}} <: AbstractMapping{Dim1,S1,S2}
    mapping::U
    field::V
    function Pairing(field::V, mapping::U) where {U,V}
        @assert codimension(mapping,2) == dimension(field) "dimensions of mapping and field are inconsistent."
        Dim1 = dimension(mapping)
        Dim2 = dimension(field)
        S1,S2 = codimension(field)
        return new{Dim1,Dim2,S1,S2,U,V}(mapping, field)
    end
end

function Base.:∘(alpha::AbstractMapping, f::GeometricMapping)
    return Pairing(alpha, f)
end

function Base.iterate(f::Pairing, state=1)
    return Pairing(f.field[state], f.mapping), state+1
end

Base.getindex(f::Pairing, k...) = Pairing(f.field[k...], f.mapping)