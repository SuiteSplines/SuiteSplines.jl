```@meta
CurrentModule = SuiteSplines
```
# SuiteSplines

[SuiteSplines](https://github.com/SuiteSplines) is a collection of Julia packages designed for applications in isogeometric analysis.

The `SuiteSplines.jl` package is intended for registration in the General repository and bundles stable and compatible versions of core packages in the SuiteSplines framework:

- IgaBase.jl
- SortedSequences.jl
- CartesianProducts.jl
- KroneckerProducts.jl
- AbstractMappings.jl
- UnivariateSplines.jl
- TensorProductBsplines.jl
- NURBS.jl
- IgaFormation.jl
- ImmersedSplines.jl
- BezierBernsteinMethods.jl
- SpecialSpaces.jl

The documentation of `SuiteSplines.jl` includes tutorials and example usage for
each of the core packages. 

## Reexporting submodules

`SuiteSplines.jl` groups core packages into submodules, such as `SuiteSplines.UnivariateSplines`.
By default, to avoid polluting the namespace, `SuiteSplines.jl` does not expose any functionality
directly. The submodules are only accessible through `SuiteSplines.UnivariateSplines`, `SuiteSplines.SortedSequences`, etc.

To conveniently access functionality from specific submodules, you can use the
[`@suitesplines_reexport`](@ref) macro.

For example, to reexport `UnivariateSplines` and `SortedSequences`, you can write:
```julia-repl
julia> using SuiteSplines

julia> @suitesplines_reexport SortedSequences UnivariateSplines

julia> I = Interval(0.0, 2.0);

julia> S = SplineSpace(2, I, 5)
SplineSpace(degree = 2, interval = [0.0, 2.0], dimension = 7)
```

If called without arguments, [`@suitesplines_reexport`](@ref) will reexport all
submodules (see [`SUITESPLINES_PKGS`](@ref)).
```julia-repl
julia> @suitesplines_reexport

julia> BezierSimplex <: BezierBernsteinMethods.AbstractSimplex
true
```
As indicated by the previous prompt, only the names originally exported by
a submodule are reexported when calling [`@suitesplines_reexport`](@ref).

## Local registry

The latest versions of core packages are registered in
[SuiteSplinesRegistry](https://github.com/SuiteSplines/SuiteSplinesRegistry).

!!! tip "Use core packages independently"
    To use the core SuiteSplines packages independently, you can add SuiteSplines registry to your depot

    ```julia-repl
    pkg> registry add https://github.com/SuiteSplines/SuiteSplinesRegistry.git
    ```


## Contributing

If you wish to contribute to SuiteSplines, please do so directly in one of the core packages. All SuiteSplines packages are listed [here](https://github.com/orgs/SuiteSplines/repositories). For more specialized, comprehensive contributions consider
contributing a new package like `TruncatedHierarchicalBsplines.jl`.