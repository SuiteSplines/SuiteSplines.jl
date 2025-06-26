using Test, SafeTestsets

@safetestset "Bezier and Bernstein methods" begin

    using BezierBernsteinMethods
    using IgaBase

    @testset "DegreeElevation" begin
        p = 3
        d = 3
        λ = BPoint(0.1,0.4,0.5)
        control = rand(dimsplinespace(p,d))
        A = degree_elevation_operator(p, d, Float64)

        @test size(A,1) == dimsplinespace(p+1,d)
        @test size(A,2) == dimsplinespace(p,d)

        # test partition of unity
        @test all(sum(A; dims=2) .≈ 1)

        # test if degree elevation preserves the polynomial
        @test bernstein(p, control, λ) ≈ bernstein(p+1, A * control, λ)
    end

    @testset "MassMatrix" begin
        p = 3
        q = p+1
        d = 4
        A = mass_matrix(Simplex(d), p, q, Float64)

        @test size(A,1) == dimsplinespace(q,d)
        @test size(A,2) == dimsplinespace(p,d)
        @test all(sum(A) ≈ 1.0/6.0) # check volume of tetrahydron

        B = mass_matrix(Simplex(d), p, p, Float64)
        x_p = rand(dimsplinespace(p,d))
        x_q = A * ((A'*A) \ (B * x_p)) # perform Galerkin projection
        @test x_q ≈ (degree_elevation_operator(p, d) * x_p)
    end
end
