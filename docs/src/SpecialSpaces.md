```@meta
CurrentModule = SuiteSplines
```
```@setup specialspaces
using SuiteSplines
@suitesplines_reexport
```
# SpecialSpaces.jl

The core function spaces used in SuiteSplines are [`SplineSpace`](@ref UnivariateSplines.SplineSpace)s and tensor-products thereof.

`SpecialSpaces.jl` introduces convenient interfaces for handling spaces used to represent scalar-, vector-valued and mixed-valued fields and geometric mappings.

## Domains and partitions

All spline spaces are defined on a [`Partition`](@ref SpecialSpaces.Partition)
of some [`Domain`](@ref SpecialSpaces.Domain).
A Cartesian product domain in two dimensions may be defined as a Cartesian product (`⨱`) of univariate [`Interval`](@ref SortedSequences.Interval)s,

```@example specialspaces
Ω = Interval(0.0, 2.0) ⨱ Interval(0.0, 4.0)
```

A partition of the domain `Ω` is just another Cartesian product of a strictly increasing points
in each dimension. To generate a uniform partition we can use [`IncreasingRange`](@ref SortedSequences.IncreasingRange),

```@example specialspaces
Δ = IncreasingRange(Ω.data[1]..., 3) ⨱ IncreasingRange(Ω.data[2], 5)
```

Less explicitly, we can also call `Partition(...)` and specify the number of
breakpoints in each dimension,

```@example specialspaces
Δ = Partition(Ω, (3,5))
```


## Scalar spline spaces

Given a partition, a tensor-product [`ScalarSplineSpace`](@ref SpecialSpaces.ScalarSplineSpace)
of degree 4 can be constructed by calling

```@example specialspaces
S = ScalarSplineSpace(2, Δ)
```

This is equivalent to
```@example specialspaces
SplineSpace(4, Δ.data[1]) ⨷ SplineSpace(2, Δ.data[2])
```

In fact, [`ScalarSplineSpace`](@ref SpecialSpaces.ScalarSplineSpace) is just a type
alias for a tensor-product of univariate spline spaces.

Conveniently, we can obtain the dimension of the tensor-product space by
```@example specialspaces
SpecialSpaces.dimension(S)
```

## Vector spline spaces

## Mixed spline spaces

### Raviart-Thomas

### Taylor-Hood

# Constraints

# Mappings