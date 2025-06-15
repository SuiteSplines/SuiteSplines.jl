# SafeTestsets does not support macros. Hence, here we
# use a module to create a safe environment
module EvaluationSetTest

using Test
using AbstractMappings, StaticArrays

@testset "EvaluationSet" begin

X = EvaluationSet{3,1}(zeros(3,4), zeros(3,4), zeros(3,4))
@test size(X) == (3,4)
@test size(X,1) == 3 && size(X,2) == 4
@test length(X) == 12
@test ndims(X) == 2
@test eltype(X) == SMatrix{3,1,Float64}
    
X = EvaluationSet{3,2}(Array{Float64}, undef, 3, 4, 5)
@test size(X) == (3,4,5)
@test size(X,1) == 3 && size(X,2) == 4 && size(X,3) == 5
@test length(X) == 60
@test ndims(X) == 3
@test eltype(X) == SMatrix{3,2,Float64}
@test elsize(X) == (3,2) && elsize(X,1)==3 && elsize(X,2)==2

Y = EvaluationSet{3,2}(zeros(3,4,5), zeros(3,4,5), zeros(3,4,5), zeros(3,4,5), zeros(3,4,5), zeros(3,4,5))
@test size(Y) == (3,4,5)
@test size(Y,1) == 3 && size(Y,2) == 4 && size(Y,3) == 5
@test length(Y) == 60
@test ndims(Y) == 3
@test eltype(Y) == SMatrix{3,2,Float64}

@test_throws AssertionError EvaluationSet{1,2}(zeros(3,4,5), zeros(3,4,2))
@test_throws AssertionError EvaluationSet{1,2}(zeros(3,4,5))

end

end