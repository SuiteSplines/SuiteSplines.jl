using Test, SafeTestsets

@safetestset "Iga formation of 2D element arrays" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U.data, V.data);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # reference values for mass and stiffness matrix
    m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v), U, V)
    m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)

    M_ref = kron(m_ref[2], m_ref[1])
    K_ref = kron(m_ref[2], m_ref_uv[1]) + kron(m_ref_uv[2], m_ref[1])

    n = prod(s -> dimsplinespace(s), U) # dimension trialspace
    m = prod(s -> dimsplinespace(s), V) # dimension testspace

    # define mass and stiffness matrix
    M = spzeros(m, n)
    K = spzeros(m, n)

    # assemble matrices
    for element in Elements(partition)

        # get global indices of test and trial functions
        A, B = TestIndices(acc, element), TrialIndices(acc, element)

        # add element mass matrix to global mass matrix
        test, trial = TestFunctions(acc, element, ders=(0,0)), TrialFunctions(acc, element, ders=(0,0))
        M[A, B] += test * trial'

        # add element stiffness matrix to global stiffness matrix
        test_u, trial_u = TestFunctions(acc, element, ders=(1,0)), TrialFunctions(acc, element, ders=(1,0))
        test_v, trial_v = TestFunctions(acc, element, ders=(0,1)), TrialFunctions(acc, element, ders=(0,1))
        K[A, B] += test_u * trial_u' + test_v * trial_v'
    end

    @testset "check_element_interval" begin
        element = Element(CartesianIndex(2,3), partition)
        dom = get_element_domain(element)
        @test dom.data[1].a == element[1,1][1] && dom.data[1].b==element[2,1][1]
        @test dom.data[2].a == element[1,1][2] && dom.data[2].b==element[1,2][2]
    end

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end
    
    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end

@safetestset "Iga formation of 2D element arrays - generalized Gaussian quadrature" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12))
    U = TensorProduct(u -> SplineSpace(3, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> GeneralizedGaussrule(u; degree=6, regularity=1), partition.data, partition.data);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # reference values for mass and stiffness matrix
    m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v), U, V)
    m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)

    M_ref = kron(m_ref[2], m_ref[1])
    K_ref = kron(m_ref[2], m_ref_uv[1]) + kron(m_ref_uv[2], m_ref[1])

    n = prod(s -> dimsplinespace(s), U) # dimension trialspace
    m = prod(s -> dimsplinespace(s), V) # dimension testspace

    # define mass and stiffness matrix
    M = spzeros(m, n)
    K = spzeros(m, n)

    # assemble matrices
    for element in Elements(partition)

        # get global indices of test and trial functions
        A, B = TestIndices(acc, element), TrialIndices(acc, element)

        # add element mass matrix to global mass matrix
        test, trial = TestFunctions(acc, element, ders=(0,0)), TrialFunctions(acc, element, ders=(0,0))
        M[A, B] += test * trial'

        # add element stiffness matrix to global stiffness matrix
        test_u, trial_u = TestFunctions(acc, element, ders=(1,0)), TrialFunctions(acc, element, ders=(1,0))
        test_v, trial_v = TestFunctions(acc, element, ders=(0,1)), TrialFunctions(acc, element, ders=(0,1))
        K[A, B] += test_u * trial_u' + test_v * trial_v'
    end

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
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # reference values for mass and stiffness matrix
    m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v), U, V)
    m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)

    M_ref = kron(m_ref[3], m_ref[2], m_ref[1])
    K_ref = kron(m_ref[3], m_ref[2], m_ref_uv[1]) + kron(m_ref[3], m_ref_uv[2], m_ref[1]) + kron(m_ref_uv[3], m_ref[2], m_ref[1])

    n = prod(s -> dimsplinespace(s), U) # dimension trialspace
    m = prod(s -> dimsplinespace(s), V) # dimension testspace

    # define mass and stiffness matrix
    M = spzeros(m, n)
    K = spzeros(m, n)

    # assemble matrices
    for element in Elements(partition)

        # get global indices of test and trial functions
        A, B = TestIndices(acc, element), TrialIndices(acc, element)

        # add element mass matrix to global mass matrix
        u = TrialFunctions(acc, element, ders=(0,0,0))
        v = TestFunctions(acc, element, ders=(0,0,0))
        M[A, B] += v * u'

        # add element stiffness matrix to global stiffness matrix
        ∇u = (TrialFunctions(acc, element, ders=(1,0,0)), TrialFunctions(acc, element, ders=(0,1,0)), TrialFunctions(acc, element, ders=(0,0,1)))
        ∇v = (TestFunctions(acc, element, ders=(1,0,0)), TestFunctions(acc, element, ders=(0,1,0)), TestFunctions(acc, element, ders=(0,0,1)))
        K[A, B] += ∇v[1] * ∇u[1]' + ∇v[2] * ∇u[2]' + ∇v[3] * ∇u[3]'
    end

    @testset "check_element_interval" begin
        element = Element(CartesianIndex(2,3,1), partition)
        dom = get_element_domain(element)
        @test dom.data[1].a == element[1,1,1][1] && dom.data[1].b==element[2,1,1][1]
        @test dom.data[2].a == element[1,1,1][2] && dom.data[2].b==element[1,2,1][2]
        @test dom.data[3].a == element[1,1,1][3] && dom.data[3].b==element[1,1,2][3]
    end

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end
    
    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end