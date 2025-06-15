using Test, SafeTestsets

@safetestset "Sorted Sequences" begin

    using SortedSequences

    @test_throws ArgumentError Interval(2,1)
    I = Interval(1.0,4.0)

    @test eltype(I) == Float64
    @test length(I) == 2
    @test size(I) == (2,)
    @test I[1] == 1.0 && I[2] == 4.0 && I[end] == 4.0
    @test_throws BoundsError I[3]

    @test measure(I) == 3.0
    @test I(2/3) ≈ 3.0

    @test_throws ArgumentError NonDecreasingVector([1,3,2,4])
    @test_throws ArgumentError IncreasingVector([1,3,2,4])
    @test_throws ArgumentError IncreasingVector([1,3,3,4])
    @test_throws ArgumentError IncreasingRange(4,1,3)
    @test_throws ArgumentError IncreasingRange(1,1,3)

    v = NonDecreasingVector([1,2,3,4])
    @test v == [1,2,3,4]
    @test length(v) == 4
    @test size(v) == (4,)

    @test v[1] == 1
    @test v[3] == 3
    @test v[4] == 4
    @test v[:] == v
    @test v[2:end-1] isa NonDecreasingVector{Int}
    @test_throws BoundsError([1, 2, 3, 4], (5,)) v[5]

    v = NonDecreasingVector([1,2,3,4])
    @test_throws ArgumentError("Values need to be in non-decreasing order.") insert!(v, 3, 1)
    @test_throws ArgumentError("Values need to be in non-decreasing order.") insert!(v, 3, 4)
    @test insert!(v, 3, 3) == [1,2,3,3,4]
    @test insert!(v, 3, 2) == [1,2,2,3,3,4]

    u = NonDecreasingVector([1.0,2.0,2.0,3.0,3.0])
    @test global_insert(u, 3) == [1.0,1.25,1.5,1.75,2.0,2.0,2.25,2.5,2.75,3.0,3.0]

    v = IncreasingVector([1,2,3,4])
    @test v == [1,2,3,4]
    @test length(v) == 4
    @test size(v) == (4,)
    
    @test IncreasingVector(v) == v
    @test all(IncreasingVector(0:4) .== [0,1,2,3,4])

    @test v[1] == 1
    @test v[3] == 3
    @test v[4] == 4
    @test v[:] == v
    @test v[2:end-1] isa IncreasingVector{Int}
    @test_throws BoundsError([1, 2, 3, 4], (5,)) v[5]

    v = IncreasingVector([1,2,4,5])
    @test_throws ArgumentError("Values need to be in increasing order.") insert!(v, 3, 2)
    @test_throws ArgumentError("Values need to be in increasing order.") insert!(v, 3, 4)
    @test insert!(v, 3, 3) == [1,2,3,4,5]

    v = NonDecreasingVector([1,1,2,3,4,4,4,5,5])
    u = Unique(v)
    @test length(u) == 5
    @test collect(u) == [(1,2), (2,1), (3,1), (4,3), (5,2)]
    u = IncreasingVector(v)
    @test u == [1,2,3,4,5]
    @test NonDecreasingVector(u, [2,1,1,3,2]) == v

    u = IncreasingVector([1.0,2.0,3.0])
    @test global_insert(u, 3) == [1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]

    u = IncreasingRange(1,5,5)
    @test u == [1,2,3,4,5]
    @test NonDecreasingVector(u, [2,1,1,3,2]) == [1,1,2,3,4,4,4,5,5]

    I = Interval(1.0,4.0)
    u = IncreasingRange(I.a, I.b, 3)
    @test length(u)==3 && u[1]==I.a && u[2]==I(0.5) && u[3]==I.b

    u = IncreasingRange(1.0,3.0,3)
    @test global_insert(u, 3) == [1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
end
