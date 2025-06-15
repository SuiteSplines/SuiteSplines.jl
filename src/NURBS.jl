module NURBS

    using LinearAlgebra, RecipesBase
    using IgaBase, AbstractMappings, CartesianProducts, KroneckerProducts
    using SortedSequences, UnivariateSplines, TensorProductBsplines

    # reexport
    export Degree, Dimension, dimension, codimension
    export Interval, IncreasingVector, IncreasingRange
    export SplineSpace, Bspline, TensorProductBspline
    export global_insert, breakpoints
    export CartesianProduct, TensorProduct, ⨱, ⨷
    export @evaluate, @evaluate!

    include("base.jl")
    include("mappings.jl")
    include("primitives.jl")
    include("refinement.jl")
    include("projection.jl")
    include("derivatives.jl")
    include("plotting.jl")

end # module
