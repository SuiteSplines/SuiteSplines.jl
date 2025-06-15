export solve!

function solve!(; displacement, mapping, distancefun, material, traction)

    # get the extended B-spline operator
    extension_operator = map(d -> spline_extension_operator(d.space, distancefun), displacement)

    # assemble system of equations
    A, b = assemble(extension_operator=extension_operator,
                displacement=displacement, 
                mapping=mapping, 
                distancefun=distancefun, 
                material=material, 
                traction=traction)

    # solve final system of equations
    u = A \ b

    # create discrete solution field
    set_solution_field!(displacement, u, extension_operator)
end

# prescribe the solution
function set_solution_field!(solution, u, E)
    dim = dimension(solution)
    s1, s2 = 0, size(E[1], 2)
    for k in 1:dim
        solution[k].coeffs[:] .= E[k] * u[s1+1:s2]
        s1 = s2
        if k<dim
            s2 += size(E[k+1], 2)
        end
    end
end