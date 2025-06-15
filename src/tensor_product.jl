export TensorProduct, ⨷

struct TensorProduct{Dim,T} <: AbstractProduct{Dim}
    data::NTuple{Dim,T}
end
TensorProduct(args...) = TensorProduct(args)
TensorProduct(f::Type, args...) = TensorProduct(map(f, args...)...)
TensorProduct(f::Function, args...) = TensorProduct(map(f, args...)...)

Base.convert(::Type{T}, x::T) where {T<:TensorProduct{1}} = x
Base.convert(::Type{<:TensorProduct{1}}, x) = TensorProduct(x)

⨷(s::T...) where T = TensorProduct(s...)
⨷(S₁::TensorProduct, s₂::T) where T = TensorProduct(S₁..., s₂)
⨷(s₁::T, S₂::TensorProduct) where T = TensorProduct(s₁, S₂...)
⨷(S₁::TensorProduct, S₂::TensorProduct) = TensorProduct(S₁..., S₂...)