export AbstractProjection, AbstractInterpolation
export project!, GalerkinProjection, Interpolation, QuasiInterpolation

abstract type AbstractProjection end
abstract type AbstractInterpolation <: AbstractProjection end

struct GalerkinProjection <:AbstractProjection end
struct Interpolation <: AbstractInterpolation end
struct QuasiInterpolation <: AbstractInterpolation end

function project! end
function project_imp! end

project!(s; onto, method::Type{<:AbstractProjection}) = map((sₖ,gₖ) -> project!(sₖ, onto=gₖ, method=method), s, onto)