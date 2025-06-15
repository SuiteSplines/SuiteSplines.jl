@recipe function f(nurbs::Nurbs{1}, ::Val{:space}; density=20)
    space = nurbs.space
    x = visualization_grid(space, density=density)
    y = bspline_interpolation_matrix(space, x, 1)[1] .* nurbs.weights'
    @evaluate w = nurbs.weightfun(x)
    return x, y./w
end