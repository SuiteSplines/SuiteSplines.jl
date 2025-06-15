module IgaFormation

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using AbstractMappings
using LinearAlgebra, LRUCache

export Interval, IncreasingRange, IncreasingVector
export CartesianProduct, TensorProduct, ⨷, ⨱
export SplineSpace, PatchRule, breakpoints, Legendre, Lobatto

include("core.jl")
include("elementform.jl")
include("elementdata.jl")
include("patchform.jl")
include("patchdata.jl")
include("sumfactcache.jl")
include("sumfactkernels.jl")

end # module
