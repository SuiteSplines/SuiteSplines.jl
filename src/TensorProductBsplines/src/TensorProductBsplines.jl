module TensorProductBsplines

    using IgaBase, SortedSequences, AbstractMappings, UnivariateSplines, CartesianProducts, KroneckerProducts

    # reexport
    export Degree, Dimension, Interval, IncreasingVector, IncreasingRange, SplineSpace
    export global_insert, breakpoints
    export CartesianProduct, TensorProduct, ⨱, ⨷
    export GeometricMapping, Field, Pairing
    export @evaluate, @evaluate!

    include("base.jl")
    include("evalcache.jl")
    include("types.jl")
    include("mappings.jl")
    include("quadrature.jl")
    include("projection.jl")
    include("derivatives.jl")
    include("refinement.jl")

end # module
