export pullback_material_law!, pullback_traction_vector!

export BMatrix, Material, material_matrix

using StaticArrays

# Hookean constitutive law
struct Material{Dim}
    youngs_modulus::Float64
    poisson_ratio::Float64
end

Material(; dim=2, youngs_modulus=1e5, poisson_ratio=0.3) = Material{dim}(youngs_modulus, poisson_ratio)

function pullback_material_law!(acc::ElementAccessor{Dim}, element::Element{Dim}, cache_nabla_data, material::Material{Dim}, mapping::GeometricMapping{Dim}) where Dim
    element_cache_nabla_data = [extract_element_cache(cache_nabla_data[i,j], element) for i in 1:Dim, j in 1:Dim]
    x = QuadraturePoints(acc, element)
    pullback_material_law!(element_cache_nabla_data, mapping, x, material)
    return element_cache_nabla_data
end

@inline function integrand_at_quad_eval(D, E_i, E_j, J)
    return (inv(J)' * E_i' * D * E_j * inv(J)) * det(J)
end

function pullback_material_law!(element_cache, mapping::GeometricMapping{Dim}, X::CartesianProduct{Dim}, material::Material{Dim}) where Dim
    
    # get the material matrix
    D = material_matrix(material)

    # get the B-matrix incidence matrices
    E = BMatrix{Dim}()

    # evaluate jacobian of the mapping at all points
    @evaluate Grad = Gradient(mapping)(X)
    
    # pull back transformation at each quadrature point
    for k in eachindex(X)

        # Jacobian at quadrature point 'k'
        J = Grad[k]

        for j in 1:Dim
            for i in 1:Dim
                element_cache[i,j][k] = integrand_at_quad_eval(D, E.data[i], E.data[j], J)
            end
        end
    end
end

function pullback_traction_vector!(acc::ElementAccessor, element::Element, cache_traction_data, mapping, traction)
    element_cache_traction_data = extract_element_cache(cache_traction_data, element)
    x = QuadraturePoints(acc, element)
    pullback_traction_vector!(squeeze(element_cache_traction_data), mapping, squeeze(x), traction)
    return element_cache_traction_data
end

@inline function traction_vector_at_quad_eval(traction, normal)
    return traction * normal
end

function pullback_traction_vector!(traction_data, mapping::GeometricMapping, X::CartesianProduct, traction)
    
    @evaluate x = mapping(X)
    @evaluate normal = Normal(mapping)(X)

    # ToDo: build sign into the evaluation of the normal
    sign = mapping.orientation

    for k in eachindex(X)
        traction_data[k] = traction_vector_at_quad_eval(traction(x[k]...), sign * normal[k])
    end
    return traction_data
end

function material_matrix(material::Material{2})
    E = material.youngs_modulus
    ν = material.poisson_ratio
    μ = E   / (2 * (1 + ν))
    λ = ν * E / ((1 + ν) * (1-2ν))
    return @SMatrix    [λ+2μ    λ       0.0;
                        λ       λ+2μ    0.0; 
                        0.0     0.0     μ];
end

function material_matrix(material::Material{3})
    E = material.youngs_modulus
    ν = material.poisson_ratio
    μ = E   / (2 * (1 + ν))
    λ = ν * E / ((1 + ν) * (1-2ν))
    return @SMatrix    [λ+2μ    λ       λ       0.0     0.0     0.0;
                        λ       λ+2μ    λ       0.0     0.0     0.0; 
                        λ       λ       λ+2μ    0.0     0.0     0.0;
                        0.0     0.0     0.0     μ       0.0     0.0;
                        0.0     0.0     0.0     0.0     μ       0.0;
                        0.0     0.0     0.0     0.0     0.0     μ];
end

struct BMatrix{Dim,F}
    data::NTuple{Dim,F}
end

function BMatrix{2}()
    E1 = @SMatrix  [1 0; 
                    0 0;
                    0 1]
    E2 = @SMatrix  [0 0; 
                    0 1;
                    1 0]
    F = typeof(E1)
    return BMatrix{2,F}((E1, E2))
end

function BMatrix{3}()
    E1 = @SMatrix  [1 0 0; 
        0 0 0;
        0 0 0;
        0 0 0;
        0 0 1;
        0 1 0]

    E2 = @SMatrix  [0 0 0; 
        0 1 0;
        0 0 0;
        0 0 1;
        0 0 0;
        1 0 0]

    E3 = @SMatrix  [0 0 0; 
        0 0 0;
        0 0 1;
        0 1 0;
        1 0 0;
        0 0 0]
    
    F = typeof(E1)
    return BMatrix{3,F}((E1, E2, E3))
end