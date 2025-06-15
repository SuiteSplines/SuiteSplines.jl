using Test, SafeTestsets

@safetestset "Iga formation of 2D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U.data, V.data);
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # reference values for mass and stiffness matrix
    m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v), U, V)
    m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)

    M_ref = kron(m_ref[2], m_ref[1])
    K_ref = kron(m_ref[2], m_ref_uv[1]) + kron(m_ref_uv[2], m_ref[1])

    n = prod(s -> dimsplinespace(s), U) # dimension trialspace
    m = prod(s -> dimsplinespace(s), V) # dimension testspace

    test, trial = TestFunctions(acc, ders=(0,0)), TrialFunctions(acc, ders=(0,0))
    M = test * trial'

    # add element stiffness matrix to global stiffness matrix
    test_u, trial_u = TestFunctions(acc, ders=(1,0)), TrialFunctions(acc, ders=(1,0))
    test_v, trial_v = TestFunctions(acc, ders=(0,1)), TrialFunctions(acc, ders=(0,1))
    K = test_u * trial_u' + test_v * trial_v'

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end
    
    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end

@safetestset "Iga formation of 3D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12), IncreasingRange(0.0,3.0,4))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U, V);
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # reference values for mass and stiffness matrix
    m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v), U, V)
    m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)

    M_ref = kron(m_ref[3], m_ref[2], m_ref[1])
    K_ref = kron(m_ref[3], m_ref[2], m_ref_uv[1]) + kron(m_ref[3], m_ref_uv[2], m_ref[1]) + kron(m_ref_uv[3], m_ref[2], m_ref[1])

    u = TrialFunctions(acc, ders=(0,0,0))
    v = TestFunctions(acc, ders=(0,0,0))
    M = v * u'

    ∇u = (TrialFunctions(acc, ders=(1,0,0)), TrialFunctions(acc, ders=(0,1,0)), TrialFunctions(acc, ders=(0,0,1)))
    ∇v = (TestFunctions(acc, ders=(1,0,0)), TestFunctions(acc, ders=(0,1,0)), TestFunctions(acc, ders=(0,0,1)))
    K = ∇v[1] * ∇u[1]' + ∇v[2] * ∇u[2]' + ∇v[3] * ∇u[3]'

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end
    
    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end