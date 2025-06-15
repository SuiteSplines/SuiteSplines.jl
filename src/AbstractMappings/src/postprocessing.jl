export l2_error,  l2_error_relative

IgaBase.standard_quadrature_rule(f, g::GeometricMapping) = IgaBase.standard_quadrature_rule(f, g[1])
IgaBase.standard_quadrature_rule(f, g::Field) = IgaBase.standard_quadrature_rule(f, g[1])


"""
    l2_error(f; to, relative=false, quadrule=standard_quadrature_rule(f,g))

Compute the (relative) L^2 error of `f` with respect to `g`. Any `f` and `g`
will do as long as `@cartesian` and `standard_quadrature_rule(f,g)` are implemented.
"""
function l2_error(f; to, quadrule=standard_quadrature_rule(to,f), relative=false)
    if relative
        return l2_error_relative_imp(f, to, quadrule)
    else
        return l2_error_imp(f, to, quadrule)
    end
end

function l2_error_imp(f, g, quadrule)

    @evaluate F = f(quadrule.x)
    @evaluate! F -= g(quadrule.x)
    
    e = zeros(elsize(F))
    for j ∈ 1:elsize(F,2)
        for i ∈ 1:elsize(F,1)
            F.data[i,j] .*= F.data[i,j]
            e[i,j] = sqrt(IgaBase.contract(F.data[i,j], quadrule.w))
        end
    end
    return e
end

function l2_error_relative_imp(f, g, quadrule)

    @evaluate F = f(quadrule.x)
    @evaluate! F -= g(quadrule.x)
    
    e = zeros(elsize(F))
    for j ∈ 1:elsize(F,2)
        for i ∈ 1:elsize(F,1)
            F.data[i,j] .*= F.data[i,j]
            e[i,j] = sqrt(IgaBase.contract(F.data[i,j], quadrule.w))
        end
    end
    @evaluate! F = g(quadrule.x)
    r = zeros(elsize(F))
    for j ∈ 1:elsize(F,2)
        for i ∈ 1:elsize(F,1)
            F.data[i,j] .*= F.data[i,j]
            r[i,j] = sqrt(IgaBase.contract(F.data[i,j], quadrule.w))
        end
    end
    return e ./ r
end