using Test, SafeTestsets

@safetestset "Simplex" begin

    using BezierBernsteinMethods

    frame = Frame(HPoint(0.0,0.0), HPoint(1.0,0.0), HPoint(1.0,1.0))
    @test frame isa Frame{3}
    simplex = Simplex(frame)
    @test dimension(simplex) == 3
    @test barycentriccoords(simplex, [0.4; 0.4]) ≈ [0.6; 0.0; 0.4]
    @test barycentriccoords(simplex, [0.3; 0.0]) ≈ [0.7; 0.3; 0.0]
    @test barycentriccoords(simplex, [1.0; 0.7]) ≈ [0.0; 0.3; 0.7]
    @test vol(simplex) ≈ 1.0/2.0

    frame = Frame(HPoint(0.0,0.0,0.0), HPoint(1.0,0.0,0.0), HPoint(1.0,1.0,0.0), HPoint(1.0,1.0,1.0))
    @test frame isa Frame{4}
    simplex = Simplex(frame)
    @test dimension(simplex) == 4
    @test barycentriccoords(simplex, [0.4; 0.0; 0.0]) ≈ [0.6; 0.4; 0.0; 0.0]
    @test vol(simplex) ≈ 1.0/6.0
end
