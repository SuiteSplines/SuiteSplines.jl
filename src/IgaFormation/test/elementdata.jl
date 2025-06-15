using Test, SafeTestsets

@safetestset "Iga formation of 2D element arrays" begin

using IgaBase, IgaFormation, UnivariateSplines
using LinearAlgebra, SparseArrays

S = SplineSpace(2, Interval(0.0,1.0), 1)
U = V = S ⨷ S
Q = TensorProduct(S -> PatchRule(breakpoints(S); npoints=Degree(S)+1, method=Legendre), U)
acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

# ground truth
M = UnivariateSplines.system_matrix(S, S)
M_u = UnivariateSplines.system_matrix(S, S, 2, 1)
M_v = UnivariateSplines.system_matrix(S, S, 1, 2)
M_uv = UnivariateSplines.system_matrix(S, S, 2, 2)

partition = acc.partition
element = Element(CartesianIndex(1,1), partition)

@testset "properties" begin
    @test size(acc) == (1,1)
    @test size(acc,1) == 1
    @test size(acc,2) == 1
    @test length(acc) == 1
    @test element.index == CartesianIndex(1,1) 
end

@testset "element mass matrix" begin
    u = TrialFunctions(acc, element; ders=(0,0))
    v = TestFunctions(acc, element; ders=(0,0)) 
    @test isapprox(norm(v * u' - kron(M, M)), 0, atol=1e-12)
end

@testset "element matrix trial_ders = (1,0)" begin    
    u = TrialFunctions(acc, element; ders=(1,0))
    v = TestFunctions(acc, element; ders=(0,0)) 
    @test isapprox(norm(v * u' - kron(M, M_u)), 0, atol=1e-12)
end

@testset "element matrix trial_ders = (0,1)" begin    
    u = TrialFunctions(acc, element; ders=(0,1))
    v = TestFunctions(acc, element; ders=(0,0)) 
    @test isapprox(norm(v * u' - kron(M_u, M)), 0, atol=1e-12)
end

@testset "element matrix test_ders = (1,0)" begin    
    u = TrialFunctions(acc, element; ders=(0,0))
    v = TestFunctions(acc, element; ders=(1,0)) 
    @test isapprox(norm(v * u' - kron(M, M_v)), 0, atol=1e-12)
end

@testset "element matrix test_ders = (0,1)" begin    
    u = TrialFunctions(acc, element; ders=(0,0))
    v = TestFunctions(acc, element; ders=(0,1)) 
    @test isapprox(norm(v * u' - kron(M_v, M)), 0, atol=1e-12)
end

@testset "element stiffness matrix term 1" begin    
    u = TrialFunctions(acc, element; ders=(1,0))
    v = TestFunctions(acc, element; ders=(1,0)) 
    @test isapprox(norm(v * u' - kron(M, M_uv)), 0, atol=1e-12)
end

@testset "element stiffness matrix term 2" begin
    u = TrialFunctions(acc, element; ders=(0,1))
    v = TestFunctions(acc, element; ders=(0,1)) 
    @test isapprox(norm(v * u' - kron(M_uv, M)), 0, atol=1e-12)
end

top    = Element(CartesianIndex(1), @view partition[1, :])
bottom = Element(CartesianIndex(1), @view partition[end, :])
left   = Element(CartesianIndex(1), @view partition[:, 1])
right  = Element(CartesianIndex(1), @view partition[:, end])

@testset "boundary elements" begin
    @test TrialIndices(acc, element) == TrialIndices(acc, top)
    @test TrialIndices(acc, element) == TrialIndices(acc, bottom)
    @test TrialIndices(acc, element) == TrialIndices(acc, left)
    @test TrialIndices(acc, element) == TrialIndices(acc, right)

    @test TestIndices(acc, element) == TestIndices(acc, top)
    @test TestIndices(acc, element) == TestIndices(acc, bottom)
    @test TestIndices(acc, element) == TestIndices(acc, left)
    @test TestIndices(acc, element) == TestIndices(acc, right)


    @test IgaFormation.get_element_indices(element) == (1,1)
    @test IgaFormation.get_element_indices(top) == (0,1)
    @test IgaFormation.get_element_indices(bottom) == (2,1)
    @test IgaFormation.get_element_indices(left) == (1,0)
    @test IgaFormation.get_element_indices(right) == (1,2)
end

end


@safetestset "Iga formation of 3D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    S = SplineSpace(2, Interval(0.0,1.0), 1)
    U = V = S ⨷ S ⨷ S
    Q = TensorProduct(S -> PatchRule(breakpoints(S); npoints=Degree(S)+1, method=Legendre), U);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # ground truth
    M = UnivariateSplines.system_matrix(S, S)
    M_u = UnivariateSplines.system_matrix(S, S, 2, 1)
    M_v = UnivariateSplines.system_matrix(S, S, 1, 2)
    M_uv = UnivariateSplines.system_matrix(S, S, 2, 2)

    partition = acc.partition
    element = Element(CartesianIndex(1,1,1), partition)
        
    @testset "properties" begin
        @test size(acc) == (1,1,1)
        @test size(acc,1) == 1
        @test size(acc,2) == 1
        @test size(acc,3) == 1
        @test length(acc) == 1
        @test element.index == CartesianIndex(1,1,1) 
    end

    @testset "element mass matrix" begin
        u = TrialFunctions(acc, element; ders=(0,0,0))
        v = TestFunctions(acc, element; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M)), 0, atol=1e-12)
    end
        
    @testset "element matrix trial_ders = (1,0,0)" begin    
        u = TrialFunctions(acc, element; ders=(1,0,0))
        v = TestFunctions(acc, element; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_u)), 0, atol=1e-12)
    end

    @testset "element matrix trial_ders = (0,1,0)" begin    
        u = TrialFunctions(acc, element; ders=(0,1,0))
        v = TestFunctions(acc, element; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M_u, M)), 0, atol=1e-12)
    end
    
    @testset "element matrix trial_ders = (0,0,1)" begin    
        u = TrialFunctions(acc, element; ders=(0,0,1))
        v = TestFunctions(acc, element; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M_u, M, M)), 0, atol=1e-12)
    end
    
    @testset "element matrix test_ders = (1,0,0)" begin    
        u = TrialFunctions(acc, element; ders=(0,0,0))
        v = TestFunctions(acc, element; ders=(1,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_v)), 0, atol=1e-12)
    end
    
    @testset "element matrix test_ders = (0,1,0)" begin    
        u = TrialFunctions(acc, element; ders=(0,0,0))
        v = TestFunctions(acc, element; ders=(0,1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_v, M)), 0, atol=1e-12)
    end
    
    @testset "element matrix test_ders = (0,0,1)" begin    
        u = TrialFunctions(acc, element; ders=(0,0,0))
        v = TestFunctions(acc, element; ders=(0,0,1)) 
        @test isapprox(norm(v * u' - kron(M_v, M, M)), 0, atol=1e-12)
    end

    @testset "element stiffness matrix term 1" begin    
        u = TrialFunctions(acc, element; ders=(1,0,0))
        v = TestFunctions(acc, element; ders=(1,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_uv)), 0, atol=1e-12)
    end
    
    @testset "element stiffness matrix term 2" begin
        u = TrialFunctions(acc, element; ders=(0,1,0))
        v = TestFunctions(acc, element; ders=(0,1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_uv, M)), 0, atol=1e-12)
    end
    
    @testset "element stiffness matrix term 2" begin
        u = TrialFunctions(acc, element; ders=(0,0,1))
        v = TestFunctions(acc, element; ders=(0,0,1)) 
        @test isapprox(norm(v * u' - kron(M_uv, M, M)), 0, atol=1e-12)
    end

    s11 = Element(CartesianIndex(1,1), @view partition[1, :, :])
    s12 = Element(CartesianIndex(1,1), @view partition[end, :, :])
    
    s21 = Element(CartesianIndex(1,1), @view partition[:, 1, :])
    s22 = Element(CartesianIndex(1,1), @view partition[:, end, :])
    
    s31 = Element(CartesianIndex(1,1), @view partition[:, :, 1])
    s32 = Element(CartesianIndex(1,1), @view partition[:, :, end])
    
    @testset "boundary elements" begin
        @test TrialIndices(acc, element) == TrialIndices(acc, s11)
        @test TrialIndices(acc, element) == TrialIndices(acc, s12)
        @test TrialIndices(acc, element) == TrialIndices(acc, s21)
        @test TrialIndices(acc, element) == TrialIndices(acc, s22)
        @test TrialIndices(acc, element) == TrialIndices(acc, s31)
        @test TrialIndices(acc, element) == TrialIndices(acc, s32)
    
        @test TestIndices(acc, element) == TestIndices(acc, s11)
        @test TestIndices(acc, element) == TestIndices(acc, s12)
        @test TestIndices(acc, element) == TestIndices(acc, s21)
        @test TestIndices(acc, element) == TestIndices(acc, s22)
        @test TestIndices(acc, element) == TestIndices(acc, s31)
        @test TestIndices(acc, element) == TestIndices(acc, s32)   
    
        @test IgaFormation.get_element_indices(element) == (1,1,1)
        @test IgaFormation.get_element_indices(s11) == (0,1,1)
        @test IgaFormation.get_element_indices(s12) == (2,1,1)
        @test IgaFormation.get_element_indices(s21) == (1,0,1)
        @test IgaFormation.get_element_indices(s22) == (1,2,1)
        @test IgaFormation.get_element_indices(s31) == (1,1,0)
        @test IgaFormation.get_element_indices(s32) == (1,1,2)
 
    end
    
end