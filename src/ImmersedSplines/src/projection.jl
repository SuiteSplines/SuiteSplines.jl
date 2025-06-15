export project!, l2_error

trialspace(b::TensorProductBspline) = b.space
trialspace(b::Field) = trialspace(b[1])

function IgaBase.project!(distancefun::AlgoimCallLevelSetFunction, field::TensorProductBspline{Dim}; onto, method::Type{GalerkinProjection}) where Dim

    # spline discretization
    U = trialspace(field)
    p = map(u -> Degree(u), U)
    partition = CartesianProduct(u -> breakpoints(u), U)

    # define regular quadrature rule, element accessor, and element sum-factorization cache
    Q = TensorProduct((d, u) -> PatchRule(d; npoints=ceil(Int, Degree(u)+1), method=Legendre), partition, U)
    acc = ElementAccessor(testspace=U, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=true);
    mat_cache = MatrixSumfactoryCache(acc);
    vec_cache = VectorSumfactoryCache(acc);

    # define cutcell quadrature rule
    Q_immersed = CutcellQuadratureRule(partition=partition, mapping=distancefun, npoints=max(p...)+1);
    acc_immersed = ElementAccessor(testspace=U, trialspace=U, quadrule=Q_immersed);

    # allocate space for global mass matrix
    m = prod(u -> dimsplinespace(u), U)
    M = spzeros(m, m)
    F = zeros(m)

    # loop over elements
    for element in Elements(partition)

        # inside-outside-test
        etest = InsideOutsideTest(distancefun, element)

        if !is_outside(etest)

            # get sumfactorization objects
            ∫ₘ = Sumfactory(mat_cache, element)  
            ∫ᵥ = Sumfactory(vec_cache, element)
                
            # case 1: element is fully inside fictituous domain
            if is_inside(etest)

                # get quadrature points
                X = QuadraturePoints(acc, element)
                Y = ones(size(X))

                # get qudrature points
                x = QuadraturePoints(acc, element)

                # get test and trial functions evaluated at quadrature points
                u = TrialFunctions(acc, element; ders=ntuple(k->0, Dim))
                v = TestFunctions(acc, element; ders=ntuple(k->0, Dim))

                # compute element mass matrix
                Mₑ = ∫ₘ(u, v, data=Y)

                # evaluate rhs forcing at quadrature nodes
                @evaluate! Y = onto(x)
                
                # compute element force vector
                fₑ = ∫ᵥ(v; data=Y)
            end

            # case 2: element intersects fictituous boundary
            if on_interface(etest)

                # get qudrature points and weights
                X = QuadraturePoints(acc_immersed, element)
                W = QuadratureWeights(acc_immersed, element)

                # loop over quadrature points
                Mₑ = ∫ₘ.data[1]; Mₑ .= 0
                fₑ = ∫ᵥ.data[1]; fₑ .= 0
                for k in eachindex(W)

                    # get quadrature point and weight
                    x, w = X[k], W[k]

                    # get test and trial functions
                    u = TrialFunctions(acc_immersed, element, k; ders=ntuple(k->0, Dim))
                    v = TestFunctions(acc_immersed, element, k; ders=ntuple(k->0, Dim))

                    # compute contribution to element mass matrix
                    Mₑ += (v .* w) * u'  

                    # compute contribution to element force vector
                    fₑ += (v .* w) * onto(x...)
                end

            end
            
            # add to global mass matrix
            A = TrialIndices(acc_immersed, element)
            B = TrialIndices(acc_immersed, element)
            M[A, B] += Mₑ
            F[A] += fₑ
        end
    end

    # compute extended splines stabilization operator
    E = spline_extension_operator(U, distancefun)

    # assmeble final system
    Mₛ = E' * M * E
    Fₛ = E' * F

    # solve system and update spline coeffs
    field.coeffs[:] = E * (Mₛ \ Fₛ) 
end

function IgaBase.standard_quadrature_rule(f, g::TensorProductBspline, distancefun::AlgoimCallLevelSetFunction)
    p = map(u -> Degree(u), g.space)
    partition = CartesianProduct(u -> breakpoints(u), U)
    return CutcellQuadratureRule(partition=partition, mapping=distancefun, npoints=max(p...)+1);
end

function AbstractMappings.l2_error(distancefun::AlgoimCallLevelSetFunction, field::AbstractMapping{Dim}; to, relative=false) where Dim

    # spline discretization
    U = trialspace(field)
    p = map(u -> Degree(u), U)
    partition = CartesianProduct(u -> breakpoints(u), U)

    # define regular quadrature rule, element accessor, and element sum-factorization cache
    Q = TensorProduct((d, u) -> PatchRule(d; npoints=ceil(Int, Degree(u)+1), method=Legendre), partition, U)
    acc = ElementAccessor(testspace=U, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=false);
   
    # define cutcell quadrature rule
    Q_immersed = CutcellQuadratureRule(partition=partition, mapping=distancefun, npoints=max(p...)+1);
    acc_immersed = ElementAccessor(testspace=U, trialspace=U, quadrule=Q_immersed);

    # loop over elements
    r = zeros(codimension(field)) # L² norm of the error
    f = zeros(codimension(field)) # L² norm of the exact solution
    for element in Elements(partition)

        # inside-outside-test
        etest = InsideOutsideTest(distancefun, element)

        # case 1: element is fully inside fictituous domain
        if !is_outside(etest)

            # get global indices
            A = TrialIndices(acc_immersed, element)

            # get sum factorization cache and input array Y
            X = QuadraturePoints(acc, element)

            # case a: element is inside
            if is_inside(etest)

                # get trial functions evaluated at quadrature points
                u = TrialFunctions(acc, element; ders=ntuple(k->0,Dim))

                # get quadrature points and weights
                x = QuadraturePoints(acc, element)
                w = QuadratureWeights(acc, element)

                # compute square error over element
                @evaluate Y = to(x)
                for k in eachindex(Y.data)
                    f[k] += contract(Y.data[k].^2, w)
                end

                # residual
                @evaluate! Y -= field(x)

                for k in eachindex(Y.data)
                    Y.data[k] .*= Y.data[k]
                    r[k] += contract(Y.data[k], w)
                end
            end

            # case 2: element intersects fictituous boundary
            if on_interface(etest)

                # get qudrature points and weights
                X = QuadraturePoints(acc_immersed, element)
                W = QuadratureWeights(acc_immersed, element)

                # loop over quadrature points
                for k in eachindex(W)

                    # compute square error over element
                    x, w = CartesianProduct(x -> [x], X[k]), W[k]
   
                    # the following is slow, but general
                    @evaluate Y = to(x)
                    for k in eachindex(Y.data)
                        f[k] += w * Y.data[k][1]^2
                    end 
   
                    # evaluate residual
                    @evaluate! Y -= field(x)

                    # l2 error at quadrature point
                    for k in eachindex(Y.data)
                        r[k] += w * Y.data[k][1]^2
                    end 
                end
            end

        end # if !is_outside(etest)

    end # for element in Elements(partition)

    return (relative) ? sqrt.(r) ./ sqrt.(f) : sqrt.(r)
end