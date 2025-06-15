using Test, SafeTestsets

@safetestset "Test ResizableArray" begin

using AbstractMappings

# create resizeable array
c = ResizableArray{Float64,3}(undef, (6,4,2))
@test size(c)==(6,4,2)
@test length(c)==48
@test eltype(c)==Float64
ptr_save = pointer(c.data)

@testset "resize without reallocation" begin    
    resize!(c, (2,3,4))
    @test size(c)==(2,3,4)
    @test length(c)==24
    @test eltype(c)==Float64
    @test pointer(c.data)==ptr_save
end

@testset "resize again without reallocation" begin    
    resize!(c, (2,12,2))
    @test size(c)==(2,12,2)
    @test length(c)==48
    @test eltype(c)==Float64
    @test pointer(c.data)==ptr_save
end

@testset "getindex and setindex!" begin
    c[1,1,1] = 1.0
    @test c[1,1,1] == 1.0
    c[1,:,end] .= 2.0
    @test all(c[1,:,end] .== 2.0)
end

end