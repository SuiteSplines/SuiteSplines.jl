using Test, SafeTestsets

@safetestset "Bezier basis evaluation" begin

    using BezierBernsteinMethods
    import BezierBernsteinMethods: step_decasteljau!
    using IgaBase

    @test grevillepoint((1,2,3)) ≈ [1.0/6.0, 2.0/6.0, 3.0/6.0]

    # test partition of unity of Bernstein polynomials
    λ = BPoint(0.1,0.4,0.5)
    let s = 0.0
        for μ in MultiIndices(3,3)
            s += bernsteinfuns(μ, λ)
        end
        @test s ≈ 1
    end

    # test single step decasteljau algorithm
    control = ones(3)
    step_decasteljau!(control, BPoint(0.1,0.4,0.5), MultiIndices(0,3))
    @test control ≈ ones(3)

    # test if decasteljau algorithm delivers same result as directly using
    # the Bernstein polynomials
    control = rand(dimsplinespace(3,3))
    b = bernstein(3, control, λ)
    @test decasteljau!(3, control, λ) ≈ b
end
