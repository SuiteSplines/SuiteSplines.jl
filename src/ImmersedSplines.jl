module ImmersedSplines

using LinearAlgebra, SparseArrays, StaticArrays
using IgaBase, AbstractMappings, SortedSequences, UnivariateSplines, KroneckerProducts, CartesianProducts
using TensorProductBsplines, IgaFormation, Algoim

include("algoim_interface.jl")
include("elementform.jl")
include("activefuns.jl")
include("projection.jl")

end # module 