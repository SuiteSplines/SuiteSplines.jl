export visualization_grid

using RecipesBase

# overload visualization_grid for your own datatypes in order to
# get a custom plotting grid. Several back-up methods have been
# implemented, but these may be inefficient for your data-type
function visualization_grid end

function AbstractMappings.visualization_grid(f::AbstractMapping{1}; density=100)
    return IncreasingRange(domain(f), density)
end

function AbstractMappings.visualization_grid(f::AbstractMapping{2}; density=(10,100))
    return CartesianProduct((x,n) -> IncreasingRange(x, n), domain(f), density)
end

@recipe function f(fun::AbstractMapping{1}; seriestype=:path)
    return fun, Val(seriestype)
end

# plot recipe spatial curves
@recipe function f(curve::AbstractMapping{1,1,1}, ::Val; density=10)
    label -->"curve",
    linewidth   --> 1
    linecolor   --> :black
    seriestype  --> :path

    x = visualization_grid(curve, density=density)
    @evaluate y = curve(x)
    return x, y
end

# plot recipe 2D spatial curves
@recipe function f(curve::AbstractMapping{1,1,2}, ::Val; density=10)
    label -->"curve",

    linewidth   --> 1
    linecolor   --> :black
    seriestype  --> :path

    x = visualization_grid(curve, density=density)
    @evaluate y = curve(x)
    return y.data[1], y.data[2]
end

# plot recipe 3D spatial curves
@recipe function f(curve::AbstractMapping{1,1,3}, ::Val; density=10)
    label -->"curve",

    linewidth   --> 1
    linecolor   --> :black
    seriestype  --> :path

    x = visualization_grid(curve, density=density)
    @evaluate y = curve(x)
    return y.data[1], y.data[2], y.data[3]
end

@recipe function f(surf::AbstractMapping{2}; seriestype=:surface)
    return surf, Val(seriestype)
end

function surface_plot_eval(f::AbstractMapping{2,1,1}, density)
    x = visualization_grid(f, density=density)
    @evaluate y = f(x)
    return x.data[2], x.data[1], y
end

function surface_plot_eval(f::AbstractMapping{2,1,2}, density)
    x = visualization_grid(f, density=density)
    @evaluate y = f(x)
    return y.data[2], y.data[1], zeros(size(y))
end

function surface_plot_eval(f::AbstractMapping{2,1,3}, density)
    x = visualization_grid(f, density=density)
    @evaluate y = f(x)
    return y.data[1], y.data[2], y.data[3] 
end

# plot recipe 2D function
@recipe function f(surface::AbstractMapping{2,1}, ::Val{S}; density=(10,10)) where S
    seriestype := S
    return surface_plot_eval(surface, density)
end

# plot recipe 2D planar surface
@recipe function f(surface::AbstractMapping{2,1,2}, ::Val{:surface}; density=(10,10))
    seriestype := :surface
    fillcolor  --> :lightgreen
    seriesalpha  --> 0.7
    legend --> :none
    return surface_plot_eval(surface, density)
end

# plot recipe 3D surface
@recipe function f(surface::AbstractMapping{2,1,3}, ::Val{:surface}; density=(10,10))
    seriestype := :surface
    fillcolor  --> :lightgreen
    seriesalpha  --> 0.7
    legend --> :none
    return surface_plot_eval(surface, density)
end

@recipe function f(mapping::AbstractMapping{3}, ::Val{:surface}; density=(10,10))
    for b in Boundary(mapping)
        @series begin
            density --> density
            b
        end
    end
end