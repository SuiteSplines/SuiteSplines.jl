export assemble

function assemble(; extension_operator, displacement::Field{Dim}, mapping::GeometricMapping{Dim}, distancefun, material::Material{Dim}, traction) where Dim

    # spline trial space
    space = displacement[1].space

    # extract partition from the solution space
    partition = CartesianProduct(s -> breakpoints(s), space)

    # define regular quadrature rule, element accessor, and element sum-factorization cache
    Q_regular = TensorProduct((d, u) -> PatchRule(d; npoints=ceil(Int, Degree(u)+1), method=Legendre), partition, space);
    acc = ElementAccessor(testspace=space, trialspace=space, quadrule=Q_regular, incorporate_weights_in_testfuns=true);

    # create cache for sum factorization
    mat_cache = MatrixSumfactoryCache(acc);
    vec_cache = VectorSumfactoryCache(acc);

    # define cutcell quadrature rule
    Q_immersed = CutcellQuadratureRule(partition=partition, mapping=distancefun, npoints=max(map(u -> Degree(u), space)...)+1);
    acc_immersed = ElementAccessor(testspace=space, trialspace=space, quadrule=Q_immersed);

    # compute the material matrix
    pullback_c!(e) = pullback_material_law!(acc, e, [FieldEvaluationCache(acc, Dim, Dim) for i in 1:Dim, j in 1:Dim], material, mapping)

    # allocate space for stiffness matrices and forcing vectors
    m = prod(s -> size(s.C,1), space)
    K = [spzeros(m, m) for i in 1:Dim, j in 1:Dim]
    f = [zeros(m) for i in 1:Dim]

    # pre-allocate element matrices
    n = prod(s -> Degree(s)+1, space)
    Kₑ = [zeros(n, n) for i in 1:Dim, j in 1:Dim]
    Fₑ = [zeros(n) for i in 1:Dim]

    # loop over elements
    for element in Elements(partition)

        # inside-outside-test
        etest = InsideOutsideTest(distancefun, element)

        if !is_outside(etest)

            # reset element arrays
            for j in 1:Dim
                for i in 1:Dim
                    Kₑ[i,j] .= 0.0
                end
                Fₑ[j] .= 0.0
            end

            # case 1: element is fully inside fictituous domain
            if is_inside(etest)
        
                # compute trial- and test-function derivatives at the quadrature points
                ∇u = map(k -> TrialFunctions(acc, element; ders=indicator(Dim, k)), 1:Dim)
                ∇v = map(k -> TestFunctions(acc, element; ders=indicator(Dim, k)), 1:Dim)

                # compute the pulled-back material data
                D_e = pullback_c!(element)

                # compute element stiffness matrices for component directions i ∈ 1:2 and j ∈ 1:2
                # these are embarrisingly parallel
                # symmetry has not been used
                for j in 1:Dim
                    for i in 1:Dim
                        # get sumfactorization object
                        ∫ = Sumfactory(mat_cache, element)
                        
                        # compute components related to first derivatives of test and trial functions
                        for beta in 1:Dim
                            for alpha in 1:Dim    
                                ∫(∇u[beta], ∇v[alpha]; data=D_e[i,j].data[alpha,beta], reset=false)
                            end
                        end

                        # add to global stiffness matrix
                        Kₑ[i,j] .= ∫.data[1]
                    
                    end  # spatial directions
                end  # spatial directions

            end # is_inside(etest)

            # case 2: element intersects fictituous boundary
            if on_interface(etest)

                # get quadrature points and weights (phase=-1 == inside)
                Q = QuadratureRule(acc_immersed, element; phase=-1)

                # get the B-matrix incidence matrices
                E = BMatrix{Dim}()
                C = material_matrix(material)

                # compute element stiffness matrices for component directions i ∈ 1:2 and j ∈ 1:2
                # these are embarrisingly parallel
                # symmetry has not been used
                for j in 1:Dim
                    for i in 1:Dim
                
                        # loop over quadrature points
                        for k in eachindex(Q.w)

                            # get quadrature point and weight
                            x, w = Q.x[k], Q.w[k]

                            # compute trial- and test-function derivatives at the quadrature points
                            ∇u = map(d -> TrialFunctions(acc_immersed, element, k; ders=indicator(Dim, d)), 1:Dim)
                            ∇v = map(d -> TestFunctions(acc_immersed, element, k; ders=indicator(Dim, d)), 1:Dim)
                            
                            # constitutive data at quadrature point
                            D_k = w * Elasticity.integrand_at_quad_eval(C, E.data[i], E.data[j], SMatrix{Dim,Dim}(I))

                            for beta in 1:Dim
                                for alpha in 1:Dim
                                    Kₑ[i,j] += (∇v[alpha] * D_k[alpha,beta])  * ∇u[beta]'
                                end
                            end
                        end # k
                    
                    end # i
                end # j

                # perform integration over boundary (phase=0)
                
                # get quadrature points and weights
                Q = QuadratureRule(acc_immersed, element; phase=0)

                # loop over quadrature points
                for k in eachindex(Q.w)

                    # get quadrature point and weight
                    x, w = Q.x[k], Q.w[k]

                    # compute trial- and test-function derivatives at the quadrature points
                    v = TestFunctions(acc_immersed, element, k; ders=ntuple(k->0, Dim))
                    
                    # constitutive data at quadrature points
                    h = traction(x...) * normal(distancefun, x)
                    for i in 1:Dim
                        Fₑ[i] += v * (w * h[i])
                    end
    
                end # k in eachindex(w)

            end # on_interface(etest)

            # add to global stiffness matrix
            A = TestIndices(acc, element)
            B = TrialIndices(acc, element)
            for j in 1:Dim
                for i in 1:Dim
                    K[i,j][A, B] += Kₑ[i,j]
                end
                f[j][A] += Fₑ[j]
            end

        end # !is_outside(etest)

    end # element loop


    # loop over non-fictitious boundaries
    for side in 1:2Dim

        # get the boundary
        mapping_∂ = boundary(mapping, side)

        # define pull-back of traction vector on current boundary component
        pullback_t!(e) = pullback_traction_vector!(acc, e, FieldEvaluationCache(acc, Dim, 1),  mapping_∂, traction)

        # pullback on immersed part takes one quadrature point at a time
        pullback_i!(x) = pullback_traction_vector!(EvaluationSet{Dim,1}(zeros(1), zeros(1), zeros(1)), mapping_∂, CartesianProduct(ξ -> [ξ], x), traction)

        # loop over boundary elements
        for element in Elements(restriction(partition, side))
        
            # inside-outside-test
            etest = InsideOutsideTest(distancefun, element)

            if !is_outside(etest)

                # reset element vector
                for i in 1:Dim
                    Fₑ[i] .= 0.0
                end

                # case 1: element is fully inside fictituous domain
                if is_inside(etest)

                    # compute test function derivatives at the quadrature points
                    v = TestFunctions(acc, element; ders=ntuple(k->0, Dim))

                    # compute the pulled-back material data
                    F_t = pullback_t!(element)

                    # loop over the spatial directions
                    for i in 1:Dim

                        #get sumfactorization object
                        ∫ = Sumfactory(vec_cache, element)

                        # compute contribution to force vector
                        ∫(v; data=F_t.data[i], reset=false)
                        
                        # add to local element vector
                        Fₑ[i] += ∫.data[1]

                    end # spatial directions
                    
                end # is_inside(etest)

                # case 2: element intersects fictituous boundary
                # 
                if on_interface(etest) && Dim>2 # critical case has not been implemented in AlgoimWraper.jl
                # for quadrature on 1D intervals

                    # get quadrature rule
                    Q = QuadratureRule(acc_immersed, element; phase=-1)
 
                    # loop over quadrature points
                    for k in eachindex(Q.w)

                        # get quadrature point and weight
                        x, w = Q.x[k], Q.w[k]

                        # compute trial- and test-function derivatives at the quadrature points
                        v = TestFunctions(acc_immersed, element, k; ders=ntuple(k->0, Dim))
                        
                        # constitutive data at quadrature points
                        h = pullback_i!(x)
                        for i in 1:Dim
                            Fₑ[i] += v * (w * h[1][i])
                        end
        
                    end # k in eachindex(w)


                end # on_interface(etest)

                # get global indices
                A = TestIndices(acc, element)

                # add to global element vector
                for i in 1:Dim
                    f[i][A] += Fₑ[i]
                end

            end # !is_outside(etest)



        end # element loop

    end # all boundaries

    # assemble global system matrix
    A, b = assemble_block_matrix_and_rhs(displacement, K, f, extension_operator)

    # return stiffness matrix and forcing vector
    return A, b
end

# indicator function
indicator(dim, k::Int) = ntuple(i -> i==k ? 1 : 0, dim)

function assemble_block_matrix_and_rhs(displacement::Field{2}, K, F, E)
    C = map(k -> sparse(KroneckerProduct(s -> s.C, displacement[k].space; reverse=true)) * E[k], 1:2)

    A = [C[1]'*K[1,1]*C[1]   C[1]'*K[1,2]*C[2]; 
         C[2]'*K[2,1]*C[1]   C[2]'*K[2,2]*C[2]]
    b = [C[1]'*F[1]; C[2]'*F[2]]
    return A, b
end

function assemble_block_matrix_and_rhs(displacement::Field{3}, K, F, E)
    C = map(k -> sparse(KroneckerProduct(s -> s.C, displacement[k].space; reverse=true)) * E[k], 1:3)

    A = [C[1]'*K[1,1]*C[1]    C[1]'*K[1,2]*C[2]     C[1]'*K[1,3]*C[3]; 
         C[2]'*K[2,1]*C[1]    C[2]'*K[2,2]*C[2]     C[2]'*K[2,3]*C[3];
         C[3]'*K[3,1]*C[1]    C[3]'*K[3,2]*C[2]     C[3]'*K[3,3]*C[3]]
    b = [C[1]'*F[1]; C[2]'*F[2]; C[3]'*F[3]]
    return A, b
end
    
ders(dim, k::Int) = ntuple(i -> i==k ? 1 : 0, dim)