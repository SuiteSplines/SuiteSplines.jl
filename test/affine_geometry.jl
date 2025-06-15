using Test, SafeTestsets

@safetestset "Affine geometry" begin

    using BezierBernsteinMethods

    @testset "Barycentric affine geometry" begin


        @test BPoint(0.2,0.3,0.5) isa BPoint{3,Float64}
        @test BVector(0.4,0.2,-0.6) isa BVector{3,Float64}
        @test_throws ArgumentError("Coordinates must sum to T(1).") BPoint(1.0,1.0,1.0)
        @test_throws ArgumentError("Coordinates must sum to T(0).") BVector(1.0,1.0,1.0)

        point = [0.2,0.3,0.5]
        λ = BPoint(point...)
        @test λ[:] == point
        @test λ[1:3] == point
        @test λ[1:2] == point[1:2]
        @test λ[1] == point[1]
        @test λ[end] == point[end]

        vec = [0.4,0.2,-0.6]
        v = BVector(vec...)
        @test v[:] == vec
        @test v[1:3] == vec
        @test v[1:2] == vec[1:2]
        @test v[1] == vec[1]
        @test v[end] == vec[end]
    end

    @safetestset "Homogenized affine geometry" begin

        using BezierBernsteinMethods

        point = [0.2,0.3,0.5]

        @test HVector(point...) isa HVector{4,Float64}
        @test HPoint(point...) isa HPoint{4,Float64}

        @test_throws  ArgumentError HVector((point...,))
        @test_throws  ArgumentError HPoint((point...,))

        p = HPoint(point...)
        @test p[:] == [point...,1.]
        @test p[1:3] == point
        @test p[1:2] == point[1:2]
        @test p[1] == point[1]
        @test p[end] == 1.0

        vec = [0.2,0.3,0.5]
        v = HVector(vec...)
        @test v[:] == [vec...,0.0]
        @test v[1:3] == vec
        @test v[1:2] == vec[1:2]
        @test v[1] == vec[1]
        @test v[end] == 0.0

        @test HPoint(0.2,0.3,0.5)-HPoint(0.2,0.3,0.3) isa HVector{4,Float64}
        @test HPoint(0.2,0.3,0.5)+HVector(0.2,0.3,0.3) isa HPoint{4,Float64}

        @test HPoint(0.2,0.3,0.5)-HPoint(0.2,0.3,0.3) ≈ [0.0,0.0,0.2,0.0]
        @test HPoint(0.2,0.3,0.5)+HVector(0.2,0.3,0.3) ≈ [0.4,0.6,0.8,1.0]
    end
end
