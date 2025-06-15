module Elasticity

using LinearAlgebra, SparseArrays, StaticArrays
using IgaBase, AbstractMappings, SortedSequences, UnivariateSplines, KroneckerProducts, CartesianProducts
using TensorProductBsplines, IgaFormation, Algoim

using ImmersedSplines

include("constitutive.jl")
include("assemble.jl")
include("solve.jl")
include("benchmarks.jl")
include("postprocessing.jl")

end # module
