module TestImmersedFormation2D

using Test

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using IgaFormation, Algoim, ImmersedSplines
using LinearAlgebra, SparseArrays, StaticArrays

# circle
phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2]]
        )

@testset "element mass matrix" begin
    
    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.0,2), IncreasingRange(0.0,1.0,2));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]);
    V = SplineSpace(4, partition.data[1]) ⨷ SplineSpace(5, partition.data[2]);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=6);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q);

    # allocate space for global mass matrix
    M = zeros(30, 12)

    # get element
    element = Element(acc,1)

    # loop over quadrature points
    qr = QuadratureRule(acc, element)

    for k in 1:length(qr)

        # get test and trial functions
        u = TrialFunctions(acc, element, k; ders=(0,0))
        v = TestFunctions(acc, element, k; ders=(0,0))

        # compute contribution to element mass matrix
        w = qr.w[k]
        M .+= (v .* w) * u'  
    end

    @test isapprox(sum(M), pi/4, atol=1e-3)
end

@testset "global mass matrix" begin
    
    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.0,10), IncreasingRange(0.0,1.0,10));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[1]);
    V = SplineSpace(4, partition.data[1]) ⨷ SplineSpace(5, partition.data[1]);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=6);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q);

    # allocate space for global mass matrix
    m = prod(v -> dimsplinespace(v), V)
    n = prod(u -> dimsplinespace(u), U)
    M = spzeros(m, n)

    # allocate element mass matrix
    Mₑ = zeros(30, 12)

    for element in Elements(partition)

        etest = InsideOutsideTest(phi, element)
        if (!is_outside(etest))

            # loop over quadrature points
            quadrule = QuadratureRule(acc, element)
            Mₑ .= 0
            for k in 1:length(quadrule)

                # get test and trial functions
                u = TrialFunctions(acc, element, k; ders=(0,0))
                v = TestFunctions(acc, element, k; ders=(0,0))

                # compute contribution to element mass matrix
                w = quadrule.w[k]
                Mₑ .+= (v .* w) * u'  
            end

            # add contribution to global mass matrix
            A = TestIndices(acc, element)
            B = TrialIndices(acc, element)
            M[A,B] += Mₑ
        end
    end

    @test isapprox(sum(M), pi/4, atol=1e-7)
end

end # module TestImmersedFormation2D


module TestImmersedFormation3D

using Test

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using IgaFormation, Algoim, ImmersedSplines
using LinearAlgebra, SparseArrays

phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0, 
        (x) -> [2.0 * x[1], 2.0 * x[2], 2.0 * x[3]]
        )

@testset "element mass matrix" begin
    
    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.0,2), IncreasingRange(0.0,1.0,2),  IncreasingRange(0.0,1.0,2));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]) ⨷ SplineSpace(4, partition.data[3]);
    V = SplineSpace(3, partition.data[1]) ⨷ SplineSpace(4, partition.data[2]) ⨷ SplineSpace(5, partition.data[3]);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=5);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q);

    # allocate space for global mass matrix
    m = prod(v -> Degree(v)+1, V)
    n = prod(u -> Degree(u)+1, U)
    M = zeros(m, n)

    # get element
    element = Element(acc,1)

    # loop over quadrature points
    qr = QuadratureRule(acc, element)
    for k in 1:length(qr)
        # get test and trial functions
        u = TrialFunctions(acc, element, k; ders=(0,0,0))
        v = TestFunctions(acc, element, k; ders=(0,0,0))

        # compute contribution to element mass matrix
        w = qr.w[k]
        M .+= (v .* w) * u'  
    end

    @test isapprox(sum(M), (4pi/3) / 8, atol=1e-3) 
end

@testset "global mass matrix" begin
        
    # spline discretization
    partition = CartesianProduct(IncreasingRange(0.0,1.0,10), IncreasingRange(0.0,1.0,11), IncreasingRange(0.0,1.0,12));
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]) ⨷ SplineSpace(4, partition.data[3]);
    V = SplineSpace(3, partition.data[1]) ⨷ SplineSpace(4, partition.data[2]) ⨷ SplineSpace(5, partition.data[3]);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=5);
    acc = ElementAccessor(testspace=V, trialspace=U, quadrule=Q);

    # allocate space for global mass matrix
    m = prod(v -> dimsplinespace(v), V)
    n = prod(u -> dimsplinespace(u), U)
    M = spzeros(m, n)

    # allocate element mass matrix
    m = prod(v -> Degree(v)+1, V)
    n = prod(u -> Degree(u)+1, U)
    Mₑ = zeros(m, n)

    for element in Elements(partition)

        etest = InsideOutsideTest(phi, element)
        if (!is_outside(etest))

            # loop over quadrature points
            quadrule = QuadratureRule(acc, element)
            Mₑ .= 0
            for k in 1:length(quadrule)

                # get test and trial functions
                u = TrialFunctions(acc, element, k; ders=(0,0,0))
                v = TestFunctions(acc, element, k; ders=(0,0,0))

                # compute contribution to element mass matrix
                w = quadrule.w[k]
                Mₑ .+= (v .* w) * u'  
            end

            # add contribution to global mass matrix
            A = TestIndices(acc, element)
            B = TrialIndices(acc, element)
            M[A,B] += Mₑ
        end
    end

    @test isapprox(sum(M), (4pi/3) / 8, atol=1e-5) 
end

end # module TestImmersedFormation3D