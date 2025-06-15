using Test, SafeTestsets

@safetestset "Multi and linear indexing" begin

    using BezierBernsteinMethods

    @test LinearIndex(2,3,2) isa LinearIndex{2,3,2}

    for d in 2:4
        for p in 1:4
            for k in 1:binomial(p+d,d)
                L = LinearIndex(p,d,k)
                @test LinearIndex(multiindex(L)) == L
            end
        end
    end

    for (i,δ) in enumerate(MultiIndices(4))
        @test δ == ntuple(j -> Int(j==i), 4)
    end

    simplex = Simplex(3)
    @test simplex isa Simplex{3}

    for (i,μ) in enumerate(MultiIndices(simplex))
        @test μ == ntuple(j -> Int(j==i), 3)
    end

    bezier = BezierSimplex(4, 3)
    @test bezier isa BezierSimplex{4,3}

    for (i,μ) in enumerate(MultiIndices(bezier))
        @test μ == multiindex(LinearIndex(4,3,i))
    end
end
