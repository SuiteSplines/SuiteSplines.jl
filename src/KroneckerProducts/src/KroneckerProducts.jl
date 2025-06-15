module KroneckerProducts

    using LinearAlgebra, TensorOperations, IgaBase

    include("base.jl")
    include("kroneckerproduct.jl")
    include("kroneckersum.jl")
    include("kronprops.jl")
    include("factorizations.jl")
    include("contract.jl")
    include("boxproduct.jl")
    include("boxprops.jl")

end #module
