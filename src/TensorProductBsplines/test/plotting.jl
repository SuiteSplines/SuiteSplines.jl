using Test, Revise

using TensorProductBsplines

using Plots

dom = Interval(0.0,π) ⨱ Interval(0.0,π)
space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), dom, (10,12));

@testset "Plot a surface function" begin
    f = GeometricMapping(dom, (x,y) -> sin(x) * sin(y))
    fₕ = GeometricMapping(TensorProductBspline, space; codimension=1)
    project!(fₕ, onto=f, method=Interpolation)

    gr()
    plot(fₕ) # plot a surface
    plot(fₕ, seriestype=:wireframe, density=(10,100)) # plot a wireframe plot
    plot(fₕ, seriestype=:contourf, density=(100,100)) # a contour plot
end

@testset "Plot a 2D planar surface" begin
    plotly()
    dom = Interval(1.0,2.0) ⨱ Interval(0.0,π)
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), dom, (10,12));

    f = GeometricMapping(dom, (r,θ) -> r*cos(θ), (r,θ) -> r*sin(θ))
    fₕ = GeometricMapping(TensorProductBspline, space; codimension=2)
    project!(fₕ, onto=f, method=Interpolation)

    plot(fₕ, fillalpha=0.7,fillcolor=:blue, density=(10,100)) # standard surface plot
    plot(fₕ, seriestype=:wireframe, density=(10,100))
end

@testset "Plot a 3D cylinder" begin
    plotly()
    dom = Interval(0.0,2*π) ⨱ Interval(0.0,π)
    space = TensorProduct((p,d,n) -> SplineSpace(Degree(p), d, n), (2,3), dom, (10,12));

    f = GeometricMapping(domain, (θ,z) -> cos(θ), (θ,z) -> sin(θ), (θ,z) -> z)
    fₕ = GeometricMapping(TensorProductBspline, space; codimension=3)
    project!(fₕ, onto=f, method=Interpolation)

    plot(fₕ, fillalpha=0.6,fillcolor=:blue, density=(100,2)) # standard surface plot
    plot(fₕ, seriestype=:wireframe, density=(100,2))
end

