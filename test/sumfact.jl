using Test, SafeTestsets

@safetestset "Global sum factorization on 2D patches" begin

    using IgaBase, IgaFormation, UnivariateSplines, KroneckerProducts
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,1.0,3), IncreasingRange(0.0,1.0,4))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    # construct the patchrule
    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U.data, V.data);

    # construct the patch accessor
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # create a sumfactorization cache to form vectors
    vec_cache = VectorSumfactoryCache(acc)
    mat_cache = MatrixSumfactoryCache(acc)

    # get the element quadrature rule
    x = QuadraturePoints(acc)

    # get the test and trial functions
    test, trial = TestFunctions(acc, ders=(0,0)), TrialFunctions(acc, ders=(0,0))

    # get evaluation cache
    cache = FieldEvaluationCache(acc, 1, 1);
    y = extract_patch_cache(cache)
    y .= 1.0

    # evaluate vector using sum factorization
    ∫ = Sumfactory(vec_cache)
    @test sum(∫(test; data=y.data[1])) ≈ 1.0
    
    # integrate the elements of the mass matrix
    ∫ = Sumfactory(mat_cache)
    M = ∫(trial, test; data=y.data[1])

    # check if the elements of the mass matrix are correctly integrated
    M_ref = KroneckerProduct((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U,V; reverse=true)
    @test isapprox(norm(M_ref-M), 0, atol=1e-12)
end

@safetestset "Global sum factorization on 3D patches" begin

    using IgaBase, IgaFormation, UnivariateSplines, KroneckerProducts
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,1.0,3), IncreasingRange(0.0,1.0,4), IncreasingRange(0.0,1.0,5))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    # construct the patchrule
    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U.data, V.data);

    # construct the patch accessor
    acc = PatchAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);

    # create a sumfactorization cache to form vectors
    vec_cache = VectorSumfactoryCache(acc)
    mat_cache = MatrixSumfactoryCache(acc)
    
    # get the element quadrature rule
    x = QuadraturePoints(acc)

    # get the test and trial functions
    test, trial = TestFunctions(acc, ders=(0,0,0)), TrialFunctions(acc, ders=(0,0,0))

    # get evaluation cache
    cache = FieldEvaluationCache(acc, 1, 1);
    y = extract_patch_cache(cache)
    y .= 1.0

    # integrate the function f(x) = 1.0 against the test functions
    ∫ = Sumfactory(vec_cache)
    @test sum(∫(test; data=y.data[1])) ≈ 1.0

    # integrate the elements of the mass matrix
    ∫ = Sumfactory(mat_cache)
    M = ∫(trial, test; data=y.data[1])

    # check if the elements of the mass matrix are correctly integrated
    M_ref = KroneckerProduct((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U,V; reverse=true)
    @test isapprox(norm(M_ref-M), 0, atol=1e-12)
end


@safetestset "Sum factorization single 2D element" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # define test and trialspace
    u1, v1 = SplineSpace(2, Interval(0.0,1.0), 1), SplineSpace(3, Interval(0.0,1.0), 1)
    u2, v2 = SplineSpace(4, Interval(0.0,2.0), 1), SplineSpace(5, Interval(0.0,2.0), 1)
    U, V = u1 ⨷ u2, v1 ⨷ v2

    # define patch quadrature rule
    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int,(Degree(u)+Degree(v))/2)+1, method=Legendre), U, V)

    # define element accessor and cache for sum factorization
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);
    cache = MatrixSumfactoryCache(acc);

    # access sumfactorization cache on element 1
    element = Element(acc, 1)
    integral = Sumfactory(cache, element)

    # compute data
    x = QuadraturePoints(acc, element)
    y = ones(size(x))

    @testset "sum-factorization mass matrix" begin
        m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U,V)
        M_ref = kron(m_ref[2], m_ref[1])

        u = TrialFunctions(acc, element; ders=(0,0))
        v = TestFunctions(acc, element; ders=(0,0))    
        M = integral(u, v; data=y)
        @test isapprox(norm(M_ref-M), 0, atol=1e-12)
    end

    @testset "sum-factorization stiffness matrix" begin
        m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U,V)
        m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U,V)
        K_ref = kron(m_ref[2], m_ref_uv[1]) + kron(m_ref_uv[2], m_ref[1])

        u = TrialFunctions(acc, element; ders=(1,0))
        v = TestFunctions(acc, element; ders=(1,0))
        integral(u, v; data=y)

        u = TrialFunctions(acc, element; ders=(0,1))
        v = TestFunctions(acc, element; ders=(0,1))

        # with reset=false, the results will be added to the previous 
        # results
        K = integral(u, v; data=y, reset=false)
        @test isapprox(norm(K_ref-K), 0, atol=1e-12)
    end

end

@safetestset "Sum factorization single 3D element" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # define test and trialspace
    U = SplineSpace(2, Interval(0.0,1.0), 1) ⨷ SplineSpace(3, Interval(0.0,2.0), 1) ⨷ SplineSpace(4, Interval(0.0,3.0), 1)
    V = SplineSpace(4, Interval(0.0,1.0), 1) ⨷ SplineSpace(5, Interval(0.0,2.0), 1) ⨷ SplineSpace(6, Interval(0.0,3.0), 1)

    # define patch quadrature rule
    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int,(Degree(u)+Degree(v))/2)+1, method=Legendre), U, V)

    # define element accessor and cache for sum factorization
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);
    cache = MatrixSumfactoryCache(acc);

    # get element 1
    element = Element(acc, 1)

    # compute data
    x = QuadraturePoints(acc, element)
    y = ones(size(x))

    @testset "sum-factorization of 3D element mass matrix" begin
        integral = Sumfactory(cache, element)
        m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U, V)
        M_ref = kron(m_ref[3], m_ref[2], m_ref[1])

        u = TrialFunctions(acc, element; ders=(0,0,0))
        v = TestFunctions(acc, element; ders=(0,0,0))    
        M = integral(u, v; data=y)
        @test isapprox(norm(M_ref-M), 0, atol=1e-12)
    end

    @testset "sum-factorization of 3D element stiffness matrix" begin
        integral = Sumfactory(cache, element)
        m_ref = map((u,v) -> UnivariateSplines.system_matrix(u, v, 1, 1), U, V)
        m_ref_uv = map((u,v) -> UnivariateSplines.system_matrix(u, v, 2, 2), U, V)
        K_ref = kron(m_ref[3], m_ref[2], m_ref_uv[1]) + kron(m_ref[3], m_ref_uv[2], m_ref[1]) + kron(m_ref_uv[3], m_ref[2], m_ref[1])

        ∇u = (TrialFunctions(acc, element; ders=(1,0,0)), TrialFunctions(acc, element; ders=(0,1,0)), TrialFunctions(acc, element; ders=(0,0,1)))
        ∇v = (TestFunctions(acc, element; ders=(1,0,0)), TestFunctions(acc, element; ders=(0,1,0)), TestFunctions(acc, element; ders=(0,0,1)))
        
        # with reset=false, the results will be added to the previous 
        # results
        local K
        for i in 1:3
            K = integral(∇u[i], ∇v[i]; data=y, reset=false)
        end

        @test isapprox(norm(K_ref - K), 0, atol=1e-12)
    end

end

@safetestset "Sum factorization of 2D mass and stiffness" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U.data, V.data);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);
    cache = MatrixSumfactoryCache(acc);

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

        integral = Sumfactory(cache, element)
        
        x = QuadraturePoints(acc, element)
        y = ones(size(x))

        # get global indices of test and trial functions
        A, B = TestIndices(acc, element), TrialIndices(acc, element)

        # add element mass matrix to global mass matrix
        test, trial = TestFunctions(acc, element, ders=(0,0)), TrialFunctions(acc, element, ders=(0,0))
        M[A, B] += integral(trial, test, data=y)

        # add element stiffness matrix to global stiffness matrix
        test_u, trial_u = TestFunctions(acc, element, ders=(1,0)), TrialFunctions(acc, element, ders=(1,0))
        test_v, trial_v = TestFunctions(acc, element, ders=(0,1)), TrialFunctions(acc, element, ders=(0,1))
        K[A, B] += integral(trial_u, test_u; data=y)
        K[A, B] += integral(trial_v, test_v; data=y)
    end

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end

    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end

@safetestset "Sum factorization of 3D mass and stiffness" begin

    using IgaBase, IgaFormation, UnivariateSplines
    using LinearAlgebra, SparseArrays

    # initialize spaces, quadrature rule and element accessor
    partition = CartesianProduct(IncreasingRange(0.0,5.0,8), IncreasingRange(0.0,6.0,12), IncreasingRange(0.0,4.0,5))
    U = TensorProduct(u -> SplineSpace(2, u), partition)
    V = TensorProduct(u -> SplineSpace(3, u), partition)

    Q = TensorProduct((u, v) -> PatchRule(breakpoints(u); npoints=ceil(Int, (Degree(u)+Degree(v))/2)+1, method=Legendre), U, V);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);
    cache = MatrixSumfactoryCache(acc);

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

        x = QuadraturePoints(acc, element)
        y = ones(size(x))

        # get global indices of test and trial functions
        A, B = TestIndices(acc, element), TrialIndices(acc, element)

        # add element mass matrix to global mass matrix
        integral = Sumfactory(cache, element)
        u = TrialFunctions(acc, element, ders=(0,0,0))
        v = TestFunctions(acc, element, ders=(0,0,0))
        M[A, B] += integral(u, v, data=y)

        # add element stiffness matrix to global stiffness matrix
        ∇u = (TrialFunctions(acc, element, ders=(1,0,0)), TrialFunctions(acc, element, ders=(0,1,0)), TrialFunctions(acc, element, ders=(0,0,1)))
        ∇v = (TestFunctions(acc, element, ders=(1,0,0)), TestFunctions(acc, element, ders=(0,1,0)), TestFunctions(acc, element, ders=(0,0,1)))

        # compute element stiffness matrix
        integral = Sumfactory(cache, element)
        for i in 1:3
            integral(∇u[i], ∇v[i]; data=y, reset=false)
        end
        K[A, B] += integral.data[1]
        
    end

    @testset "global mass matrix" begin
        @test isapprox(norm(M - M_ref), 0.0, atol=1e-12)
    end

    @testset "global stiffness matrix" begin
        @test isapprox(norm(K - K_ref), 0.0, atol=1e-12)
    end

end


@safetestset "Sum factorization low level" begin

    using IgaFormation

    @testset "3-d case - element array" begin
        m, n = (3,4,5), (4,5,6)
        dims = (7,7,7)
        trialfuns = ntuple(k -> ones(Float64, n[k], dims[k]), length(n))
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)


        A, B = pre_allocate_arrays(Float64, dims, m, n)
        sumfact!(A, B, C, testfuns, trialfuns)
        @test all(A[1] .== sum(C))
    end

    @testset "3d case - element vector" begin
        m, n = (3,4,5), (1,1,1)
        dims = (7,7,7)
        trialfuns = ntuple(k -> ones(Float64, dims[k])', length(n))
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)

        A, B = pre_allocate_arrays(Float64, dims, m, n)
        sumfact!(A, B, C, testfuns, trialfuns)
        @test all(A[1] .== sum(C))
    end

    @testset "3d case - rhs vector" begin
        m  = (3,4,5)
        dims = (7,7,7)
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)

        A = pre_allocate_arrays(Float64, dims, m)
        sumfact!(A, C, testfuns)
        @test all(A[1] .== sum(C))
    end

    @testset "2-d case - element array" begin
        m, n = (3,4), (4,5)
        dims = (7,7)
        trialfuns = ntuple(k -> ones(Float64, n[k], dims[k]), length(n))
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)

        A, B = pre_allocate_arrays(Float64, dims, m, n)
        sumfact!(A, B, C, testfuns, trialfuns)
        @test all(A[1] .== sum(C))
    end

    @testset "2-d case - element vector" begin
        m, n = (3,4), (1,1)
        dims = (7,7)
        trialfuns = ntuple(k -> ones(Float64, dims[k])', length(n))
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)

        A, B = pre_allocate_arrays(Float64, dims, m, n)
        sumfact!(A, B, C, testfuns, trialfuns)
        @test all(A[1] .== sum(C))
    end

    @testset "2-d case - rhs vector" begin
        m  = (3,4)
        dims = (7,7)
        testfuns = ntuple(k -> ones(Float64, m[k], dims[k]), length(m))
        C = ones(dims...)

        A = pre_allocate_arrays(Float64, dims, m)
        sumfact!(A, C, testfuns)
        @test all(A[1] .== sum(C))
    end
    
    @testset "3-d case - general stiffness matrix" begin
        q = (10,10,10)
        c = ones(q...)
        
        m = (3,3,3)
        n = (3,3,3)
        
        tfuns = map((mₖ, qₖ) -> rand(mₖ, qₖ), m, q)
        
        A, B = pre_allocate_arrays(Float64, q, m, n)
        sumfact!(A, B, c, tfuns, tfuns)
        
        @test isapprox(A[1], transpose(A[1]), atol=eps(Float64))
        
        Aref = (c[:] .* kron(tfuns[3]', kron(tfuns[2]', tfuns[1]')))' * kron(tfuns[3]',kron(tfuns[2]', tfuns[1]'))
        @test isapprox(Aref, A[1], atol=1e-12)
    end

    @testset "3-d case - general rhs vector" begin
        q = (10,10,10)
        c = ones(q...)
        
        m = (3,3,3)
        
        tfuns = map((mₖ, qₖ) -> rand(mₖ, qₖ), m, q)
        
        A = pre_allocate_arrays(Float64, q, m)
        sumfact!(A, c, tfuns)
        
        Aref = sum(kron(tfuns[3]', kron(tfuns[2]', tfuns[1]')) .* c[:], dims=1)'
        @test isapprox(Aref, A[1], atol=1e-12)
    end

    @testset "2-d case - general stiffness matrix" begin
        q = (10,10)
        c = rand(q...)
        
        m = (3,3)
        n = (3,3)
        
        tfuns = map((mₖ, qₖ) -> rand(mₖ, qₖ), m, q)
        
        A, B = pre_allocate_arrays(Float64, q, m, n)
        sumfact!(A, B, c, tfuns, tfuns)
        
        @test isapprox(A[1], transpose(A[1]), atol=eps(Float64))
        
        Aref = (c[:] .* kron(tfuns[2]', tfuns[1]'))' * kron(tfuns[2]', tfuns[1]')
        @test isapprox(Aref, A[1], atol=1e-12)
    end

    @testset "2-d case - general rhs vector" begin
        q = (10,10)
        c = ones(q...)
        
        m = (3,3)
        
        tfuns = map((mₖ, qₖ) -> rand(mₖ, qₖ), m, q)
        
        A = pre_allocate_arrays(Float64, q, m)
        sumfact!(A, c, tfuns)
        
        Aref = sum(kron(tfuns[2]', tfuns[1]') .* c[:], dims=1)'
        @test isapprox(Aref, A[1], atol=1e-12)
    end

end