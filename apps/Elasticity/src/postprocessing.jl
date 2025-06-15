export CauchyStress

struct CauchyStress{Dim,T1<:GeometricMapping{Dim}, T2<:AbstractMapping{Dim}} <: AbstractMapping{Dim,Dim,Dim}
    mapping::T1
    deformation::T2
    material::Material{Dim}
    function CauchyStress(mapping::GeometricMapping{Dim}, deformation::AbstractMapping{Dim}, material::Material{Dim}) where Dim
        T1 = typeof(mapping)
        T2 = typeof(deformation)
        new{Dim,T1,T2}(mapping, deformation, material)
    end
end

IgaBase.standard_quadrature_rule(f, g::CauchyStress) = IgaBase.standard_quadrature_rule(f, g.deformation)

ImmersedSplines.trialspace(s::CauchyStress) = ImmersedSplines.trialspace(s.deformation)

# enables fast evaluation of stress
for op in [:(=), :(-=)]

    local S = Val{op}

    for Dim in 2:3

        @eval function AbstractMappings.evalkernel!(op::$S, y::EvaluationSet{$Dim,$Dim}, x::CartesianProduct{$Dim}, σ::CauchyStress{$Dim})

            # initialize
            mapping = σ.mapping
            displacement = σ.deformation
            material = σ.material  

            # get the material properties
            E = material.youngs_modulus
            ν = material.poisson_ratio
            μ = E / (2.0 * (1.0 + ν))
            λ = ν * E / ((1.0 + ν) * (1.0-2.0ν))
            
            # identity matrix
            Id = SMatrix{$Dim,$Dim}(I)

            # evaluate jacobian of the mapping at all points
            # and the gradient of displacement
            @evaluate ∂F = Gradient(mapping)(x)
            @evaluate ∂u = Gradient(displacement)(x)
            
            # pull back transformation at each quadrature point
            for k in eachindex(x)
                update_at_evaluation_point!(op, y, ∂F[k], ∂u[k], Id, λ, μ, k)        
            end
            return y
        end

    end # Dim
end # OP

@inline function update_at_evaluation_point!(::Val{:(=)}, y, ∂F, ∂u, Id, λ, μ, k)
    ∇u = inv(∂F) * ∂u
    ε = 0.5 * (∇u + ∇u')
    y[k] = λ * tr(ε) * Id + 2.0 * μ * ε        
end

@inline function update_at_evaluation_point!(::Val{:(-=)}, y, ∂F, ∂u, Id, λ, μ, k)
    ∇u = inv(∂F) * ∂u
    ε = 0.5 * (∇u + ∇u')
    y[k] -= λ * tr(ε) * Id + 2.0 * μ * ε        
end