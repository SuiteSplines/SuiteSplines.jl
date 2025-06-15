using Test, SafeTestsets

@safetestset "Array utilities" begin

    using SortedSequences

    import SortedSequences.count_multiplicity_down
    import SortedSequences.count_multiplicity_up

    @test construct_vector([0.0,1.0,3.0], [3,2,3]) == [0.0,0.0,0.0,1.0,1.0,3.0,3.0,3.0]

    # count the number of times that v[index] is repeated, starting
    # from index
    v = [1,1,2,3,4,4,4,5,5]
    @test count_multiplicity_down(v, 1) == 2
    @test count_multiplicity_down(v, 2) == 1
    @test count_multiplicity_down(v, 3) == 1
    @test count_multiplicity_down(v, 4) == 1
    @test count_multiplicity_down(v, 5) == 3
    @test count_multiplicity_down(v, 6) == 2
    @test count_multiplicity_down(v, 7) == 1
    @test count_multiplicity_down(v, 8) == 2
    @test count_multiplicity_down(v, 9) == 1

    @test count_multiplicity_up(v, 1) == 1
    @test count_multiplicity_up(v, 2) == 2
    @test count_multiplicity_up(v, 3) == 1
    @test count_multiplicity_up(v, 4) == 1
    @test count_multiplicity_up(v, 5) == 1
    @test count_multiplicity_up(v, 6) == 2
    @test count_multiplicity_up(v, 7) == 3
    @test count_multiplicity_up(v, 8) == 1
    @test count_multiplicity_up(v, 9) == 2

    v = [0,0,1,2,2,3,3,3,4,5,5]
    u, m, ia, ic = deconstruct_vector(v)
    @test length(u) == 6
    @test u == [0,1,2,3,4,5]
    @test m == [2,1,2,3,1,2]
    @test v == construct_vector(u,m)
    @test u[ia] == v
    @test v[ic] == u
end
