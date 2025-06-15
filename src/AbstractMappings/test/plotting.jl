using Test, Revise

using SortedSequences, CartesianProducts, AbstractMappings

using Plots

@testset "Plot a regular function" begin
    domain = Interval(0.0,3.0)
    f = GeometricMapping(domain, x -> sin(x) * x^2 + x + 1)

    plot(f; linewidth=2, density=100)
end

@testset "Plot a 3/4 circle" begin
    domain = Interval(0.0,1.5*π)
    f = GeometricMapping(domain, x -> 2 * cos(x), x -> 2 * sin(x))

    plot(f; linewidth=2, density=100)
end

@testset "Plot a 3D helix" begin
    domain = Interval(0.0,8*π)
    f = GeometricMapping(domain, x -> 2 * cos(x), x -> 2 * sin(x), x -> x/(2*π))

    plot(f; linecolor=:red, linewidth=2, density=200)
end

@testset "Plot a surface function" begin
    domain = Interval(0.0,π) ⨱ Interval(0.0,π)
    f = GeometricMapping(domain, (x,y) -> sin(x) * sin(y))

    plot(f) # plot a surface
    plot(f, seriestype=:wireframe, density=(10,100)) # plot a wireframe plot
    plot(f, seriestype=:contourf, density=(100,100)) # a contour plot
end

@testset "Plot a 2D planar surface" begin
    pyplot()
    domain = Interval(1.0,2.0) ⨱ Interval(0.0,π)
    f = GeometricMapping(domain, (r,θ) -> r*cos(θ), (r,θ) -> r*sin(θ))

    plot(f, fillalpha=0.1,fillcolor=:blue, density=(10,100)) # standard surface plot
    plot(f, seriestype=:wireframe, density=(10,100))
end

@testset "Plot a 3D cylinder" begin
    pyplot()
    domain = Interval(0.0,2*π) ⨱ Interval(0.0,π)
    f = GeometricMapping(domain, (θ,z) -> cos(θ), (θ,z) -> sin(θ), (θ,z) -> z)

    plot(f, fillalpha=0.1,fillcolor=:blue, density=(100,2)) # standard surface plot
    plot(f, seriestype=:wireframe, density=(100,2))
end

