module AbstractMappings

using StaticArrays, LRUCache, ForwardDiff, LinearAlgebra
using IgaBase, SortedSequences, CartesianProducts

export CartesianProduct, ⨱, IncreasingRange, IncreasingVector, Interval

import CartesianProducts: CartesianProduct, ⨱
import SortedSequences: Interval, IncreasingRange, IncreasingVector

include("array.jl")
include("evalset.jl")
include("mappings.jl")
include("fields.jl")
include("evaluate.jl")
include("plotting.jl")
include("postprocessing.jl")


end # module
