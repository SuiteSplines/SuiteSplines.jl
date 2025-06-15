module IgaBase

using LinearAlgebra, SparseArrays

include("base.jl")
include("mappings.jl")
include("projection.jl")
include("refinement.jl")
include("linalg.jl")
include("quadrature.jl")

end # module
