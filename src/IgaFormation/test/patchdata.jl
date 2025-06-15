using Test, SafeTestsets

@safetestset "Iga formation of 2D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    S = SplineSpace(2, Interval(0.0,1.0), 1)
    U = V = S ⨷ S
    Q = TensorProduct(S -> PatchRule(breakpoints(S); npoints=Degree(S)+1, method=Legendre), U)
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # ground truth
    M = UnivariateSplines.system_matrix(S, S)
    M_u = UnivariateSplines.system_matrix(S, S, 2, 1)
    M_v = UnivariateSplines.system_matrix(S, S, 1, 2)
    M_uv = UnivariateSplines.system_matrix(S, S, 2, 2)

    @testset "mass matrix" begin
        u = TrialFunctions(acc; ders=(0,0))
        v = TestFunctions(acc; ders=(0,0)) 
        @test isapprox(norm(v * u' - kron(M, M)), 0, atol=1e-12)
    end

    @testset "matrix trial_ders = (1,0)" begin    
        u = TrialFunctions(acc; ders=(1,0))
        v = TestFunctions(acc; ders=(0,0)) 
        @test isapprox(norm(v * u' - kron(M, M_u)), 0, atol=1e-12)
    end

    @testset "matrix trial_ders = (0,1)" begin    
        u = TrialFunctions(acc; ders=(0,1))
        v = TestFunctions(acc; ders=(0,0)) 
        @test isapprox(norm(v * u' - kron(M_u, M)), 0, atol=1e-12)
    end

    @testset "matrix test_ders = (1,0)" begin    
        u = TrialFunctions(acc; ders=(0,0))
        v = TestFunctions(acc; ders=(1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_v)), 0, atol=1e-12)
    end

    @testset "matrix test_ders = (0,1)" begin    
        u = TrialFunctions(acc; ders=(0,0))
        v = TestFunctions(acc; ders=(0,1)) 
        @test isapprox(norm(v * u' - kron(M_v, M)), 0, atol=1e-12)
    end

    @testset "stiffness matrix term 1" begin    
        u = TrialFunctions(acc; ders=(1,0))
        v = TestFunctions(acc; ders=(1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_uv)), 0, atol=1e-12)
    end

    @testset "stiffness matrix term 2" begin
        u = TrialFunctions(acc; ders=(0,1))
        v = TestFunctions(acc; ders=(0,1)) 
        @test isapprox(norm(v * u' - kron(M_uv, M)), 0, atol=1e-12)
    end

end

@safetestset "Iga formation of 3D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    S = SplineSpace(2, Interval(0.0,1.0), 1)
    U = V = S ⨷ S ⨷ S
    Q = TensorProduct(S -> PatchRule(breakpoints(S); npoints=Degree(S)+1, method=Legendre), U);
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # ground truth
    M = UnivariateSplines.system_matrix(S, S)
    M_u = UnivariateSplines.system_matrix(S, S, 2, 1)
    M_v = UnivariateSplines.system_matrix(S, S, 1, 2)
    M_uv = UnivariateSplines.system_matrix(S, S, 2, 2)

    partition = acc.partition
    
    @testset "mass matrix" begin
        u = TrialFunctions(acc; ders=(0,0,0))
        v = TestFunctions(acc; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M)), 0, atol=1e-12)
    end
        
    @testset "matrix trial_ders = (1,0,0)" begin    
        u = TrialFunctions(acc; ders=(1,0,0))
        v = TestFunctions(acc; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_u)), 0, atol=1e-12)
    end

    @testset "matrix trial_ders = (0,1,0)" begin    
        u = TrialFunctions(acc; ders=(0,1,0))
        v = TestFunctions(acc; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M_u, M)), 0, atol=1e-12)
    end
    
    @testset "matrix trial_ders = (0,0,1)" begin    
        u = TrialFunctions(acc; ders=(0,0,1))
        v = TestFunctions(acc; ders=(0,0,0)) 
        @test isapprox(norm(v * u' - kron(M_u, M, M)), 0, atol=1e-12)
    end
    
    @testset "matrix test_ders = (1,0,0)" begin    
        u = TrialFunctions(acc; ders=(0,0,0))
        v = TestFunctions(acc; ders=(1,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_v)), 0, atol=1e-12)
    end
    
    @testset "matrix test_ders = (0,1,0)" begin    
        u = TrialFunctions(acc; ders=(0,0,0))
        v = TestFunctions(acc; ders=(0,1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_v, M)), 0, atol=1e-12)
    end
    
    @testset "matrix test_ders = (0,0,1)" begin    
        u = TrialFunctions(acc; ders=(0,0,0))
        v = TestFunctions(acc; ders=(0,0,1)) 
        @test isapprox(norm(v * u' - kron(M_v, M, M)), 0, atol=1e-12)
    end

    @testset "stiffness matrix term 1" begin    
        u = TrialFunctions(acc; ders=(1,0,0))
        v = TestFunctions(acc; ders=(1,0,0)) 
        @test isapprox(norm(v * u' - kron(M, M, M_uv)), 0, atol=1e-12)
    end
    
    @testset "stiffness matrix term 2" begin
        u = TrialFunctions(acc; ders=(0,1,0))
        v = TestFunctions(acc; ders=(0,1,0)) 
        @test isapprox(norm(v * u' - kron(M, M_uv, M)), 0, atol=1e-12)
    end
    
    @testset "stiffness matrix term 2" begin
        u = TrialFunctions(acc; ders=(0,0,1))
        v = TestFunctions(acc; ders=(0,0,1)) 
        @test isapprox(norm(v * u' - kron(M_uv, M, M)), 0, atol=1e-12)
    end
    
end