using Test, Revise

using SortedSequences, UnivariateSplines, TensorProductBsplines, NURBS
using Plots

@testset "Project and plot a nurbs function and its space" begin
    gr() # use 'gr' as a plotting backend

    # Get the spline interpolation
    space = SplineSpace(Degree(2), Interval(0.0,pi), 10)
    nurbs = Nurbs(space)
    nurbs.weights .= 0.8.*rand(dimsplinespace(space)) .+ 0.2

    g = ScalarFunction(x -> sin.(x.^2 + x))
    project!(nurbs, onto=g, method=Interpolation)

    h1 = plot(nurbs; density=100, label="Nurbs interpolant of g(x)")
    plot!(nurbs; seriestype=:controlpolygon)
    plot!(nurbs; seriestype=:knots)
    plot!(nurbs; seriestype=:grevillepoints)

    x = IncreasingRange(0.0,pi,200)
    plot!(x, g(x); label="g(x)")

    h2 = plot(nurbs; seriestype=:space)
    plot(h1, h2, layout=(2,1))
end

@testset "Plot a Nurbs arc" begin
    gr() # use 'gr' as a plotting backend
    curve = arc(origin=(1.0,1.0), radius=1.0, α=0.0, β=2π/3)
    h1 = plot(curve; density=30,
        aspect_ratio=:equal,
        axis=:off,
        label="arc")
    plot!(curve; seriestype=:controlpolygon)
    plot!(curve; seriestype=:knots)
    h2 = plot(curve[1]; seriestype=:space)
    plot(h1, h2, layout=(2,1))
end

@testset "Plot a Nurbs circle" begin
    gr() # use 'gr' as a plotting backend
    curve = circle(origin=(1.0,1/2), radius=2)
    h1 = plot(curve; density=100, label="circle",
        aspect_ratio=:equal,
        axis=:off,
        ticks=false)
    plot!(curve; seriestype=:controlpolygon)
    plot!(curve; seriestype=:knots)
    h2 = plot(curve[1]; seriestype=:space)
    plot(h1, h2, layout=(2,1))
end

@testset "Plot a Nurbs surface" begin
    plotly() # use 'plotly' as a plotting backend for surfaces and volumes
    s = cylinder(radius=1.0, height=2.0)
    h = plot(s; density=(100,10), label="cylinder", axis=:off)
end

#@testset "Plot a Nurbs surface" begin
plotly() # use 'plotly' as a plotting backend for surfaces and volumes
s = partial_tube(inner_radius=0.25, outer_radius=0.75, height=1.0, α = 0.0, β = π)
#h = plot(size=(1200,1200));
plot(s; density=(20,20), label="cylinder", axis=:off)

boundary(s,1)


#end