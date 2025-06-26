module BezierBernsteinMethods

    using LinearAlgebra, SparseArrays, StaticArrays, Combinatorics
    using IgaBase

    include("bezier_types.jl")
    include("affine_geometry.jl")
    include("simplex.jl")
    include("multi_indexing.jl")
    include("bezierfuns.jl")
    include("bezier_methods.jl")
    include("whithney_forms.jl")

end # module
