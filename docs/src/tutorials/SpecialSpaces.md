```@meta
CurrentModule = SuiteSplines
```
```@setup specialspaces
using SuiteSplines
@suitesplines_reexport
```
```@setup specialspaces_higherdims
using SuiteSplines
@suitesplines_reexport
```
!!! warning "Experimental interface"

    This package implements interfaces that may eventually be integrated into the `IgaBase.jl` package.

# SpecialSpaces.jl

The core function spaces in SuiteSplines are [`SplineSpace`](@ref UnivariateSplines.SplineSpace) and tensor-products thereof.

`SpecialSpaces.jl` introduces convenient interfaces for handling spaces used to represent scalar-, vector- and mixed-valued fields and geometric mappings.

!!! details "Note on higher dimensions"

    Almost everything described on this page generalizes to higher dimensions, e.g.
    ```@repl specialspaces_higherdims
    Ω = Interval(0.0, 2.0) ⨱ Interval(0.0, 4.0) ⨱ Interval(-1.0, 1.0)
    Δ = Partition(Ω, (3, 5, 7));
    S = ScalarSplineSpace(4, Δ)
    ```
    The only exception are the Raviart-Thomas and Taylor-Hood spaces, which are
    implemented only in two and three dimensions.

## Domains and partitions

All spline spaces are defined on a [`Partition`](@ref SpecialSpaces.Partition)
of some [`Domain`](@ref SpecialSpaces.Domain).
A Cartesian product domain in two dimensions may be defined as a Cartesian product of intervals,

```@repl specialspaces
Ω = Interval(0.0, 2.0) ⨱ Interval(0.0, 4.0)
```

A partition of the domain `Ω` is just another Cartesian product of a strictly increasing real numbers
on the intervals defining each domain dimension. To generate an uniform partition we can use
[`IncreasingRange`](@ref SortedSequences.IncreasingRange),

```@repl specialspaces
Δ = IncreasingRange(Ω.data[1]..., 3) ⨱ IncreasingRange(Ω.data[2], 5)
```

Less explicitly, we can also call

```@repl specialspaces
Δ = Partition(Ω, (3, 5))
```
which gives us a Cartesian product partition with 3 breakpoints in the first and 5 breakpoints the second dimension.

## Scalar spline spaces

Given a partition, a tensor-product [`ScalarSplineSpace`](@ref SpecialSpaces.ScalarSplineSpace)
of degree 2 can be constructed by calling

```@repl specialspaces
S = ScalarSplineSpace(2, Δ)
```

This is equivalent to
```@repl specialspaces
SplineSpace(2, Δ.data[1]) ⨷ SplineSpace(2, Δ.data[2])
```

In fact, [`ScalarSplineSpace`](@ref SpecialSpaces.ScalarSplineSpace) is just a type
alias for a tensor-product of univariate spline spaces. The dimension of the tensor-product
was inferred from the dimension of the provided partition. `SpecialSpaces.jl`
implements a few more handy [constructors](@ref SpecialSpaces.ScalarSplineSpace)
for tensor-product spline spaces, e.g. with different degrees in each dimension.

Conveniently, we can obtain the dimension and linear indices of the tensor-product space by
```@repl specialspaces
linds = indices(S)
dim = dimension(S)
```
Linear indices are useful when handling solution vectors, fields and geometric mappings
on more complex spaces. These spaces follow next.

## Vector spline spaces

A tensor-product [`VectorSplineSpace`](@ref SpecialSpaces.VectorSplineSpace)
can be constructed from multiple scalar spline spaces,
```@repl specialspaces
V = VectorSplineSpace(S, S, S)
```
The requirement here is that the scalar spaces are defined on the same partition.

To access the scalar space corresponding to the `k`th component of a vector space, we can
index into the collection,
```@repl specialspaces
V[1]
V[1] == V[2] == V[3]
```

Note that in this example the domain dimension is two and the codomain dimension is three.
More commonly, the dimension and codimension coincide and a vector space of degree 2 can be
defined directly as
```@repl specialspaces
V = VectorSplineSpace(2, Δ)
```

The linear indices for each component are
```@repl specialspaces
for k in 1:2
    @show indices(V, k)
end
```

Note here the subtlety of calling `indices(V, k)` instead of `indices(V[k])`.
The element type of `V` is
```@repl specialspaces
eltype(V)
```
Thus an expression like the one below will call the `indices` method for a scalar spline space,
```@repl specialspaces
indices.(V)
```
This does not give the desired result. The following expression on the other hand does,
```@repl specialspaces
indices.(Ref(V), 1:2)
```
As expected, the complete range of linear indices for a vector space is
```@repl specialspaces
indices(V)
```


## Mixed spline spaces
Mixed spline function spaces are defined as custom structs subtyping
[`MixedSplineSpace`](@ref SpecialSpaces.MixedSplineSpace), where each struct
field corresponds to a scalar or vector spline space.

A basic definition of a custom mixed spline space might look like the following:
```julia
struct TaylorHood{Dim,T} <: MixedSplineSpace{Dim,T}
    V::VectorSplineSpace{Dim,Dim,T}
    Q::ScalarSplineSpace{Dim,T}
    
    function TaylorHood(p::Degree, Δ::Partition{Dim,T}) where {Dim,T<:Real}
        @assert p ≥ 2
        p = ntuple(i -> p, Dim)
        V = VectorSplineSpace(ntuple(i -> ScalarSplineSpace(p, Δ), Dim)...)
        Q = ScalarSplineSpace(p .- 1, Δ)
        new{Dim,T}(V, Q)
    end
end
```
`SpecialSpaces.jl` implements interfaces that allow us to operate on such spaces.

### Taylor–Hood
One implemented example of a mixed spline space is the
[`TaylorHood`](@ref SpecialSpaces.TaylorHood) space, which defines an inf-sub stable
pair of spaces for velocities and pressure.

```@repl specialspaces
th = TaylorHood(2, Δ)
dimension(th)
dimension(th, :Q)
dimension(th, :V)
dimension(th, :V, 1)
dimension(th, :V, 2)
indices(th, :Q)
indices(th, :V)
indices(th, :V, 1)
indices(th, :V, 2)
```

### Raviart–Thomas
Another implemented example of a mixed spline space is the
[`RaviartThomas`](@ref SpecialSpaces.RaviartThomas) space, which constructs
a structure preserving (divergence conforming) pair of spline spaces for
velocities and pressure.
```@repl specialspaces
rt = RaviartThomas(2, Δ)
dimension(rt)
dimension(rt, :Q)
dimension(rt, :V)
dimension(rt, :V, 1)
dimension(rt, :V, 2)
indices(rt, :Q)
indices(rt, :V)
indices(rt, :V, 1)
indices(rt, :V, 2)
```

Both examples are implemented for dimensions two and three and arbitrary degrees.

### Iterable mixed spline space

If a custom struct defining a mixed spline space is not desired, `SpecialSpaces.jl`
implements an iterable mixed spline space for general use.
[`IterableMixedSplineSpace`](@ref SpecialSpaces.IterableMixedSplineSpace)
can be constructed from a named tuple of scalar and vector spline spaces.
```@repl specialspaces
V = VectorSplineSpace(2, Δ)
Q = ScalarSplineSpace(1, Δ)
ith = IterableMixedSplineSpace((V=V, Q=Q))
dimension(ith, :Q)
dimension(ith, :V)
```
The example above reproduces the Taylor–Hood mixed spline space.

In addition to the usual interface for mixed spaces, it is also possible to iterate over
the spaces collected in the iterable mixed spline space,
```@repl specialspaces
map(dimension, ith)
```

## Mappings
`SpecialSpaces.jl` implements convenient constructors for fields and geometric mappings
defined on special spaces introduced above. For example,

```@repl specialspaces
ϕʰ = GeometricMapping(Ω, S);
t = Field(V);
uʰ = GeometricMapping(Ω, rt, :V);
pʰ = GeometricMapping(Ω, rt, :Q);
```

Now, given a solution vector `x` from some computation with the Raviart–Thomas
space, we can set the coefficients of the geometric mappings `uʰ` and `pʰ` as follows

```@repl specialspaces
x = rand(dimension(rt));
setcoeffs!(uʰ, rt, :V, x);
setcoeffs!(pʰ, rt, :Q, x);
```

We can also obtain vertically concatenated coefficients of fields and geometric mappings,
```@repl specialspaces
getcoeffs(uʰ) == x[indices(rt, :V)]
getcoeffs(pʰ) == x[indices(rt, :Q)] 
```

The interface for scalar and vector spaces is identical but without the field symbol.

## Constraints

Recall that space constraints on univariate [`SplineSpace`](@ref UnivariateSplines.SplineSpace)
are enforced by extraction operators,
```@repl specialspaces
space = SplineSpace(4, IncreasingRange(0.0, π, 5); cleft=[1,2])
space.C # extraction operator
```
Naturally, space constraints on tensor-product spaces can be enforced by
Kronecker product extraction operators. To avoid tedious work on the level of
univariate spline spaces, i.e. defining univariate spline spaces with
the arguments `cleft`, `cright`, `cperiodic`, `SpecialSpaces.jl` implements
containers for collecting constraints on scalar, vector and mixed spline spaces.
Upon passing these collections to respective space constructors,
the contraints are distributed down to the univariate spaces.

### Scalar spline space constraints

On the lowest abstraction level, there are three rather self-explanatory methods:

- [`left_constraint!`](@ref SpecialSpaces.left_constraint!)
- [`right_constraint!`](@ref SpecialSpaces.right_constraint!)
- [`periodic_constraint!`](@ref SpecialSpaces.periodic_constraint!)

These methods operate on scalar constraint containers and push constraint
conditions for a univariate spline space in dimension `dim` to the container.

The following is equivalent to setting `cleft=[1]` argument on
the spline space in the first dimension of the tensor-product space.
```@repl specialspaces
C = ScalarSplineSpaceConstraints{2}()
left_constraint!(C; dim=1) # or left_constraint!(S; c=[1,], dim=1)
```

This constraint container with dimension two can be then applied while
constructing a tensor-product space of dimension two. Let us reuse the
existing scalar spline space `S` and construct a scalar spline space with
constraint conditions defined in `C`
```@repl specialspaces
S
S[1].C
S = ScalarSplineSpace(S, C)
S[1].C
```

### Vector spline space constraints

Existing spaces can be used to construct suitable (empty) constraint containers.
For example, to construct a constraints container for the vector spline space `V`
we can write
```@repl specialspaces
C = VectorSplineSpaceConstraints(V)
```

[`VectorSplineSpaceConstraints`](@ref SpecialSpaces.VectorSplineSpaceConstraints)
is just a collection of
[`ScalarSplineSpaceConstraints`](@ref SpecialSpaces.ScalarSplineSpaceConstraints)
so that `left_constraint!`, `right_constraint!` and `periodic_constraint!` can be
applied to each element of this container.

```@repl specialspaces
left_constraint!(C[1]; dim=1)
V[1]
V = VectorSplineSpace(V, C);
V[1]
```

### Mixed spline space constraints

[`MixedSplineSpaceConstraints`](@ref SpecialSpaces.MixedSplineSpaceConstraints)
is just a type alias for a named tuple.

```@repl specialspaces
C = MixedSplineSpaceConstraints((V=VectorSplineSpaceConstraints(rt.V), Q=ScalarSplineSpaceConstraints(rt.Q)))
C.V
C.Q
```

### Clamped constraints

For convenience, `SuiteSplines.jl` provides
[`clamped_constraint!`](@ref SpecialSpaces.clamped_constraint!) for scalar and
vector spline spaces which can be set on edges/faces of the domain. For example,
```@repl specialspaces
C = VectorSplineSpaceConstraints(V)
clamped_constraint!(C, :top)
```
This is equivalent to setting `left_constraint!` with the arguments `(c=[1], dim=2)` in
each dimension of the vector space. It is also possible to clamp multiple edges
in a subset of dimensions,
```@repl specialspaces
C = VectorSplineSpaceConstraints{3}()
clamped_constraint!(C, :left, :right, :bottom, :top, :back, :front; dim=[1,3])
```
In this case, the scalar spline spaces in dimensions 1 and 3 are clamped on the whole
boundary of the domain. No constraints are applied in dimension 2.

Clamping for scalar spline spaces works the same way. Valid boundary labels are:
`:left`, `:right`, `:bottom`, `:top`, `:back`, `:front`, where the last two are only
valid in three dimensions.

!!! tip "Numbering and labeling boundaries on Cartesian grids"

    SuiteSplines follows a simple convention for numbering boundaries of Cartesian
    product domains: boundaries are numbered dimension by dimension by assigning
    ever increasing integers first to the beginning then to the end of the interval
    in a particular dimension. The labels `:left → :right`, `:bottom → :top`, `:back → :front`
    are directed always in the direction of the axes. See
    [`boundary_number`](@ref SpecialSpaces.boundary_number) for more details.

    

### Extraction operators

To obtain the Kronecker product extraction operator of a scalar spline space use
```@repl specialspaces
extraction_operator(S)
extraction_operator(S; sparse=true)
```

The parameter `sparse` can be used if we wish to obtain the Kronecker product
extraction operator as a sparse matrix. 

To obtain the Kronecker product extraction operator of a vector spline space use
```@repl specialspaces
extraction_operator.(V)
extraction_operator.(V; sparse=true)
```
or
```@repl specialspaces
extraction_operators(V)
extraction_operators(V; sparse=true)
```

Note, that the extraction operator is *not* stored in the constraints containers!


## API
```@autodocs
Modules = [
    SuiteSplines.SpecialSpaces,
]
Order   = [:function, :type, :macro, :constant]
```