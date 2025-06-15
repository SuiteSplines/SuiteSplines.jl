```@meta
CurrentModule = SuiteSplines
```
# SuiteSplines



## Reexporting submodules

`SuiteSplines` groups packages into submodules, such as `SuiteSplines.UnivariateSplines`.
By default, to avoid polluting the namespace, `using SuiteSplines` or `import SuiteSplines` does
not expose any functionality directly. The submodules are only accessible through
`SuiteSplines.UnivariateSplines`.

To conveniently access functionality from specific submodules, you can use the
[`@suitesplines_reexport`](@ref) macro. For example, to reexport `UnivariateSplines`
and `SortedSequences`, you can write:
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



## Index
```@index
```

```@autodocs
Modules = [
    SuiteSplines,
    SuiteSplines.IgaBase,
    SuiteSplines.SortedSequences,
    SuiteSplines.CartesianProducts,
    SuiteSplines.KroneckerProducts,
    SuiteSplines.AbstractMappings,
    SuiteSplines.BezierBernsteinMethods,
    SuiteSplines.UnivariateSplines,
    SuiteSplines.TensorProductBsplines,
    SuiteSplines.NURBS,
    SuiteSplines.IgaFormation,
    SuiteSplines.ImmersedSplines,
]
Order   = [:function, :type, :macro, :constant]
```